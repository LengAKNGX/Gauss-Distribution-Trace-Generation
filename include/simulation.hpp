#pragma once

#include <cstddef>
#include <iosfwd>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "curve.hpp"
#include "vec2.hpp"

namespace gauss_trace {

constexpr double pi = 3.141592653589793238462643383279502884;

double wrap_angle(double angle);
Vec2 heading_vector(double theta);

struct State {
    Vec2 P;
    double theta = pi / 2.0;
    double v = 1.0;
    double t = 0.0;
};

struct SimulationConfig {
    int n_steps = 1000;
    double dt = 0.05;
    double lookahead_L = 1.0;
    double alpha = 0.75;
    double beta = 0.2;
    double v_exp = 1.0;
    double sigma_theta = 0.03;
    double sigma_v = 0.02;
    double sigma_measure = 0.05;
    int sample_interval_N = 10;
    unsigned long long seed = 20260611ULL;
    double local_projection_delta = 0.25;
    int max_newton_iters = 12;
    double newton_tol = 1e-10;
    bool speed_nonnegative = true;
    ArcLengthMethod arc_length_method = ArcLengthMethod::SecondOrder;
    std::string output_csv_path = "outputs/trace.csv";
};

struct StepRecord {
    int index = 0;
    double time = 0.0;
    Vec2 true_position;
    Vec2 scan_position;
    double theta = 0.0;
    double theta_exp = 0.0;
    double v = 0.0;
    double t_min = 0.0;
    double lookahead_T = 0.0;
    double epsilon_theta = 0.0;
    double epsilon_v = 0.0;
    double measure_dx = 0.0;
    double measure_dy = 0.0;
    bool was_sampled = false;
};

double squared_distance_to_curve(const Curve2D& curve, const Vec2& point, double t);

double local_project_parameter(const Curve2D& curve,
                               const Vec2& previous_position,
                               const Vec2& next_position,
                               double previous_t,
                               double delta,
                               int max_newton_iters,
                               double newton_tol);

class Simulation {
public:
    Simulation(const Curve2D& curve, SimulationConfig config, State initial_state);

    StepRecord step(int index);
    std::vector<StepRecord> run();

    const State& state() const { return state_; }
    const SimulationConfig& config() const { return config_; }

private:
    const Curve2D& curve_;
    SimulationConfig config_;
    State state_;
    std::mt19937_64 rng_;
    std::normal_distribution<double> angle_noise_;
    std::normal_distribution<double> speed_noise_;
    std::normal_distribution<double> measure_noise_;
};

void write_csv(std::ostream& output,
               const Curve2D& curve,
               const SimulationConfig& config,
               const State& initial_state,
               const std::vector<StepRecord>& records);

void write_csv_file(const std::string& path,
                    const Curve2D& curve,
                    const SimulationConfig& config,
                    const State& initial_state,
                    const std::vector<StepRecord>& records);

}  // namespace gauss_trace
