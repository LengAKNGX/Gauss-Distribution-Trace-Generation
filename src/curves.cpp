#include "curve.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace gauss_trace {
namespace {

double simpson_integral(const Curve2D& curve, double t, double increment, int intervals) {
    if (intervals % 2 != 0) {
        ++intervals;
    }
    const double h = increment / static_cast<double>(intervals);
    double sum = curve.speed(t) + curve.speed(t + increment);
    for (int i = 1; i < intervals; ++i) {
        const double u = t + h * static_cast<double>(i);
        sum += (i % 2 == 0 ? 2.0 : 4.0) * curve.speed(u);
    }
    return sum * h / 3.0;
}

}  // namespace

ArcLengthMethod parse_arc_length_method(const std::string& value) {
    if (value == "exact") {
        return ArcLengthMethod::Exact;
    }
    if (value == "first" || value == "first_order") {
        return ArcLengthMethod::FirstOrder;
    }
    if (value == "second" || value == "second_order") {
        return ArcLengthMethod::SecondOrder;
    }
    if (value == "numerical") {
        return ArcLengthMethod::Numerical;
    }
    throw std::invalid_argument("unknown arc length method: " + value);
}

std::string to_string(ArcLengthMethod method) {
    switch (method) {
    case ArcLengthMethod::Exact:
        return "exact";
    case ArcLengthMethod::FirstOrder:
        return "first_order";
    case ArcLengthMethod::SecondOrder:
        return "second_order";
    case ArcLengthMethod::Numerical:
        return "numerical";
    }
    return "unknown";
}

double Curve2D::speed(double t) const {
    return norm(derivative(t));
}

double Curve2D::lookahead_increment(double t, double length, ArcLengthMethod method) const {
    (void)method;
    const double current_speed = speed(t);
    if (current_speed <= 0.0) {
        throw std::runtime_error("curve speed must be positive for lookahead");
    }
    return length / current_speed;
}

std::string VerticalLineCurve::name() const {
    return "vertical_line";
}

Vec2 VerticalLineCurve::position(double t) const {
    return {0.0, t};
}

Vec2 VerticalLineCurve::derivative(double) const {
    return {0.0, 1.0};
}

Vec2 VerticalLineCurve::second_derivative(double) const {
    return {0.0, 0.0};
}

double VerticalLineCurve::speed(double) const {
    return 1.0;
}

double VerticalLineCurve::lookahead_increment(double, double length, ArcLengthMethod) const {
    return length;
}

EllipseCurve::EllipseCurve(double a_value, double b_value) : a_(a_value), b_(b_value) {
    if (a_ <= 0.0 || b_ <= 0.0) {
        throw std::invalid_argument("ellipse axes must be positive");
    }
}

std::string EllipseCurve::name() const {
    return "ellipse";
}

Vec2 EllipseCurve::position(double t) const {
    return {a_ * std::cos(t), b_ * std::sin(t)};
}

Vec2 EllipseCurve::derivative(double t) const {
    return {-a_ * std::sin(t), b_ * std::cos(t)};
}

Vec2 EllipseCurve::second_derivative(double t) const {
    return {-a_ * std::cos(t), -b_ * std::sin(t)};
}

double EllipseCurve::speed(double t) const {
    const double sin_t = std::sin(t);
    const double cos_t = std::cos(t);
    return std::sqrt(a_ * a_ * sin_t * sin_t + b_ * b_ * cos_t * cos_t);
}

double EllipseCurve::lookahead_increment(double t, double length, ArcLengthMethod method) const {
    const double f = speed(t);
    if (f <= 0.0) {
        throw std::runtime_error("ellipse speed must be positive");
    }

    switch (method) {
    case ArcLengthMethod::Exact:
    case ArcLengthMethod::Numerical:
        return numerical_lookahead(t, length);
    case ArcLengthMethod::FirstOrder:
        return length / f;
    case ArcLengthMethod::SecondOrder: {
        // Second-order Lagrange-inversion approximation for ellipse arc length.
        const double delta = a_ * a_ - b_ * b_;
        return length / f - (delta * std::sin(2.0 * t) / (4.0 * std::pow(f, 4.0))) * length * length;
    }
    }
    return length / f;
}

double EllipseCurve::numerical_lookahead(double t, double length) const {
    if (length == 0.0) {
        return 0.0;
    }
    const double sign = length > 0.0 ? 1.0 : -1.0;
    const double target = std::abs(length);
    double upper = sign * std::max(target / speed(t), 1e-6);
    while (std::abs(arc_length_from(t, upper)) < target) {
        upper *= 2.0;
    }

    double low = 0.0;
    double high = upper;
    for (int i = 0; i < 80; ++i) {
        const double mid = 0.5 * (low + high);
        const double arc = std::abs(arc_length_from(t, mid));
        if (arc < target) {
            low = mid;
        } else {
            high = mid;
        }
    }
    return 0.5 * (low + high);
}

double EllipseCurve::arc_length_from(double t, double increment) const {
    return simpson_integral(*this, t, increment, 96);
}

std::unique_ptr<Curve2D> make_curve(const std::string& curve_type, double ellipse_a, double ellipse_b) {
    if (curve_type == "line" || curve_type == "vertical_line") {
        return std::make_unique<VerticalLineCurve>();
    }
    if (curve_type == "ellipse") {
        return std::make_unique<EllipseCurve>(ellipse_a, ellipse_b);
    }
    throw std::invalid_argument("unknown curve type: " + curve_type);
}

}  // namespace gauss_trace
