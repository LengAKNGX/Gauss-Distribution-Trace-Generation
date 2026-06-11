#include "simulation.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace gauss_trace {
namespace {

double stationarity(const Curve2D& curve, const Vec2& point, double t) {
    return dot(point - curve.position(t), curve.derivative(t));
}

double stationarity_derivative(const Curve2D& curve, const Vec2& point, double t) {
    const Vec2 d1 = curve.derivative(t);
    const Vec2 d2 = curve.second_derivative(t);
    return -norm_squared(d1) + dot(point - curve.position(t), d2);
}

double choose_better_parameter(const Curve2D& curve, const Vec2& point, double lhs, double rhs) {
    const double lhs_distance = squared_distance_to_curve(curve, point, lhs);
    const double rhs_distance = squared_distance_to_curve(curve, point, rhs);
    return lhs_distance <= rhs_distance ? lhs : rhs;
}

double distribution_stddev(double sigma) {
    return sigma > 0.0 ? sigma : 1.0;
}

}  // namespace

double wrap_angle(double angle) {
    double wrapped = std::fmod(angle + pi, 2.0 * pi);
    if (wrapped <= 0.0) {
        wrapped += 2.0 * pi;
    }
    return wrapped - pi;
}

Vec2 heading_vector(double theta) {
    return {std::cos(theta), std::sin(theta)};
}

double squared_distance_to_curve(const Curve2D& curve, const Vec2& point, double t) {
    return norm_squared(point - curve.position(t));
}

double local_project_parameter(const Curve2D& curve,
                               const Vec2& previous_position,
                               const Vec2& next_position,
                               double previous_t,
                               double delta,
                               int max_newton_iters,
                               double newton_tol) {
    const Vec2 tangent = curve.derivative(previous_t);
    const double tangent_length_sq = norm_squared(tangent);
    if (tangent_length_sq <= 0.0) {
        throw std::runtime_error("curve derivative vanished during projection");
    }

    const double predicted_t = previous_t + dot(next_position - previous_position, tangent) / tangent_length_sq;
    const double left = predicted_t - delta;
    const double right = predicted_t + delta;

    double best_t = choose_better_parameter(curve, next_position, left, right);

    bool converged = false;
    double t = predicted_t;
    for (int iter = 0; iter < max_newton_iters; ++iter) {
        const double denominator = stationarity_derivative(curve, next_position, t);
        if (std::abs(denominator) < 1e-14) {
            break;
        }

        const double next_t = t - stationarity(curve, next_position, t) / denominator;
        if (next_t < left || next_t > right) {
            converged = false;
            break;
        }
        if (std::abs(next_t - t) < newton_tol) {
            t = next_t;
            converged = true;
            break;
        }
        t = next_t;
    }

    if (converged && t >= left && t <= right) {
        const double newton_distance = squared_distance_to_curve(curve, next_position, t);
        const double best_distance = squared_distance_to_curve(curve, next_position, best_t);
        if (newton_distance < best_distance) {
            best_t = t;
        }
    }

    return best_t;
}

Simulation::Simulation(const Curve2D& curve, SimulationConfig config, State initial_state)
    : curve_(curve),
      config_(std::move(config)),
      state_(initial_state),
      rng_(config_.seed),
      angle_noise_(0.0, distribution_stddev(config_.sigma_theta)),
      speed_noise_(0.0, distribution_stddev(config_.sigma_v)),
      measure_noise_(0.0, distribution_stddev(config_.sigma_measure)) {
    if (config_.n_steps < 0) {
        throw std::invalid_argument("n_steps must be nonnegative");
    }
    if (config_.dt <= 0.0) {
        throw std::invalid_argument("dt must be positive");
    }
    if (config_.sample_interval_N <= 0) {
        throw std::invalid_argument("sample_interval_N must be positive");
    }
    if (config_.local_projection_delta < 0.0) {
        throw std::invalid_argument("local_projection_delta must be nonnegative");
    }
}

StepRecord Simulation::step(int index) {
    const State previous = state_;

    const double lookahead_T = curve_.lookahead_increment(previous.t, config_.lookahead_L, config_.arc_length_method);
    const Vec2 q_plus = curve_.position(previous.t + lookahead_T);

    Vec2 expected_direction;
    try {
        expected_direction = normalized(q_plus - previous.P);
    } catch (const std::runtime_error&) {
        expected_direction = normalized(curve_.derivative(previous.t));
    }
    const double theta_exp = std::atan2(expected_direction.y, expected_direction.x);

    const double epsilon_theta = config_.sigma_theta > 0.0 ? angle_noise_(rng_) : 0.0;
    const double epsilon_v = config_.sigma_v > 0.0 ? speed_noise_(rng_) : 0.0;

    state_.theta = wrap_angle(previous.theta + config_.alpha * wrap_angle(theta_exp - previous.theta) + epsilon_theta);
    state_.v = previous.v + config_.beta * (config_.v_exp - previous.v) + epsilon_v;
    if (config_.speed_nonnegative) {
        state_.v = std::max(0.0, state_.v);
    }
    state_.P = previous.P + state_.v * config_.dt * heading_vector(state_.theta);
    state_.t = local_project_parameter(curve_,
                                       previous.P,
                                       state_.P,
                                       previous.t,
                                       config_.local_projection_delta,
                                       config_.max_newton_iters,
                                       config_.newton_tol);

    const bool was_sampled = (index % config_.sample_interval_N) == 0;
    double measure_dx = 0.0;
    double measure_dy = 0.0;
    Vec2 scan_position{std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::quiet_NaN()};
    if (was_sampled) {
        measure_dx = config_.sigma_measure > 0.0 ? measure_noise_(rng_) : 0.0;
        measure_dy = config_.sigma_measure > 0.0 ? measure_noise_(rng_) : 0.0;
        scan_position = state_.P + Vec2{measure_dx, measure_dy};
    }

    return StepRecord{index,
                      static_cast<double>(index) * config_.dt,
                      state_.P,
                      scan_position,
                      state_.theta,
                      theta_exp,
                      state_.v,
                      state_.t,
                      lookahead_T,
                      epsilon_theta,
                      epsilon_v,
                      measure_dx,
                      measure_dy,
                      was_sampled};
}

std::vector<StepRecord> Simulation::run() {
    std::vector<StepRecord> records;
    records.reserve(static_cast<std::size_t>(config_.n_steps));
    for (int index = 1; index <= config_.n_steps; ++index) {
        records.push_back(step(index));
    }
    return records;
}

void write_csv(std::ostream& output,
               const Curve2D& curve,
               const SimulationConfig& config,
               const State& initial_state,
               const std::vector<StepRecord>& records) {
    output << std::setprecision(17);
    output << "# curve_type," << curve.name() << '\n';
    output << "# n_steps," << config.n_steps << '\n';
    output << "# dt," << config.dt << '\n';
    output << "# lookahead_L," << config.lookahead_L << '\n';
    output << "# alpha," << config.alpha << '\n';
    output << "# beta," << config.beta << '\n';
    output << "# v_exp," << config.v_exp << '\n';
    output << "# sigma_theta," << config.sigma_theta << '\n';
    output << "# sigma_v," << config.sigma_v << '\n';
    output << "# sigma_measure," << config.sigma_measure << '\n';
    output << "# sample_interval_N," << config.sample_interval_N << '\n';
    output << "# seed," << config.seed << '\n';
    output << "# local_projection_delta," << config.local_projection_delta << '\n';
    output << "# max_newton_iters," << config.max_newton_iters << '\n';
    output << "# newton_tol," << config.newton_tol << '\n';
    output << "# speed_nonnegative," << (config.speed_nonnegative ? "true" : "false") << '\n';
    output << "# arc_length_method," << to_string(config.arc_length_method) << '\n';
    output << "# initial_x," << initial_state.P.x << '\n';
    output << "# initial_y," << initial_state.P.y << '\n';
    output << "# initial_theta," << initial_state.theta << '\n';
    output << "# initial_v," << initial_state.v << '\n';
    output << "# initial_t," << initial_state.t << '\n';
    output << "index,time,true_x,true_y,scan_x,scan_y,theta,theta_exp,v,t_min,T,"
              "epsilon_theta,epsilon_v,measure_dx,measure_dy,was_sampled\n";

    for (const StepRecord& record : records) {
        output << record.index << ','
               << record.time << ','
               << record.true_position.x << ','
               << record.true_position.y << ','
               << record.scan_position.x << ','
               << record.scan_position.y << ','
               << record.theta << ','
               << record.theta_exp << ','
               << record.v << ','
               << record.t_min << ','
               << record.lookahead_T << ','
               << record.epsilon_theta << ','
               << record.epsilon_v << ','
               << record.measure_dx << ','
               << record.measure_dy << ','
               << (record.was_sampled ? 1 : 0) << '\n';
    }
}

void write_csv_file(const std::string& path,
                    const Curve2D& curve,
                    const SimulationConfig& config,
                    const State& initial_state,
                    const std::vector<StepRecord>& records) {
    const std::filesystem::path output_path(path);
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("failed to open output CSV: " + path);
    }
    write_csv(output, curve, config, initial_state, records);
}

}  // namespace gauss_trace
