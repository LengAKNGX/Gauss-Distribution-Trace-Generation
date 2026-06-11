#include "curve.hpp"
#include "simulation.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gauss_trace;

namespace {

bool close(double lhs, double rhs, double tolerance = 1e-10) {
    return std::abs(lhs - rhs) <= tolerance;
}

void test_wrap_angle() {
    const std::vector<double> inputs{-10.0 * pi, -pi, -0.1, 0.0, pi, 7.5};
    for (double value : inputs) {
        const double wrapped = wrap_angle(value);
        assert(wrapped > -pi);
        assert(wrapped <= pi);
    }
    assert(close(wrap_angle(-pi), pi));
}

void test_line_atan2_direction() {
    VerticalLineCurve curve;
    const Vec2 p{0.0, 4.0};
    const double lookahead = curve.lookahead_increment(4.0, 2.0, ArcLengthMethod::Exact);
    const Vec2 q_plus = curve.position(4.0 + lookahead);
    const Vec2 direction = normalized(q_plus - p);
    const double theta_exp = std::atan2(direction.y, direction.x);
    assert(close(theta_exp, pi / 2.0));
}

void test_line_one_step_moves_toward_axis() {
    VerticalLineCurve curve;
    SimulationConfig config;
    config.n_steps = 1;
    config.dt = 0.1;
    config.lookahead_L = 1.0;
    config.alpha = 1.0;
    config.beta = 1.0;
    config.v_exp = 1.0;
    config.sigma_theta = 0.0;
    config.sigma_v = 0.0;
    config.sigma_measure = 0.0;
    config.sample_interval_N = 1;
    config.arc_length_method = ArcLengthMethod::Exact;

    State right{{2.0, 0.0}, pi / 2.0, 1.0, 0.0};
    Simulation sim_right(curve, config, right);
    const StepRecord right_step = sim_right.step(1);
    assert(right_step.true_position.x < right.P.x);

    State left{{-2.0, 0.0}, pi / 2.0, 1.0, 0.0};
    Simulation sim_left(curve, config, left);
    const StepRecord left_step = sim_left.step(1);
    assert(left_step.true_position.x > left.P.x);
}

void test_reproducible_measurement_noise() {
    VerticalLineCurve curve;
    SimulationConfig config;
    config.n_steps = 5;
    config.sigma_theta = 0.01;
    config.sigma_v = 0.01;
    config.sigma_measure = 0.1;
    config.sample_interval_N = 1;
    config.seed = 12345;
    State initial{{1.0, 0.0}, pi / 2.0, 1.0, 0.0};

    Simulation lhs(curve, config, initial);
    Simulation rhs(curve, config, initial);
    const auto lhs_records = lhs.run();
    const auto rhs_records = rhs.run();
    assert(lhs_records.size() == rhs_records.size());
    for (std::size_t i = 0; i < lhs_records.size(); ++i) {
        assert(lhs_records[i].true_position.x == rhs_records[i].true_position.x);
        assert(lhs_records[i].scan_position.x == rhs_records[i].scan_position.x);
        assert(lhs_records[i].measure_dx == rhs_records[i].measure_dx);
        assert(lhs_records[i].measure_dy == rhs_records[i].measure_dy);
    }
}

void test_local_projection_endpoint_fallback() {
    EllipseCurve curve(3.0, 1.5);
    const Vec2 previous = curve.position(0.25);
    const Vec2 next{5.0, 0.2};
    const double previous_t = 0.25;
    const double projected = local_project_parameter(curve, previous, next, previous_t, 0.05, 0, 1e-12);

    const Vec2 tangent = curve.derivative(previous_t);
    const double predicted = previous_t + dot(next - previous, tangent) / norm_squared(tangent);
    const double left = predicted - 0.05;
    const double right = predicted + 0.05;
    assert(close(projected, left) || close(projected, right));
}

void test_ellipse_speed_formula() {
    EllipseCurve curve(4.0, 2.0);
    const double t = 0.7;
    const double expected = std::sqrt(4.0 * 4.0 * std::sin(t) * std::sin(t)
                                    + 2.0 * 2.0 * std::cos(t) * std::cos(t));
    assert(close(curve.speed(t), expected));
}

void test_ellipse_second_order_lookahead() {
    const double a = 5.0;
    const double b = 2.0;
    const double t = 0.9;
    const double length = 0.3;
    EllipseCurve curve(a, b);
    const double f = curve.speed(t);
    const double delta = a * a - b * b;
    const double expected = length / f - (delta * std::sin(2.0 * t) / (4.0 * std::pow(f, 4.0))) * length * length;
    const double actual = curve.lookahead_increment(t, length, ArcLengthMethod::SecondOrder);
    assert(close(actual, expected));
}

}  // namespace

int main() {
    test_wrap_angle();
    test_line_atan2_direction();
    test_line_one_step_moves_toward_axis();
    test_reproducible_measurement_noise();
    test_local_projection_endpoint_fallback();
    test_ellipse_speed_formula();
    test_ellipse_second_order_lookahead();
    std::cout << "All smoke tests passed.\n";
    return 0;
}
