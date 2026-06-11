#pragma once

#include <memory>
#include <string>

#include "vec2.hpp"

namespace gauss_trace {

enum class ArcLengthMethod {
    Exact,
    FirstOrder,
    SecondOrder,
    Numerical
};

ArcLengthMethod parse_arc_length_method(const std::string& value);
std::string to_string(ArcLengthMethod method);

class Curve2D {
public:
    virtual ~Curve2D() = default;

    virtual std::string name() const = 0;
    virtual Vec2 position(double t) const = 0;
    virtual Vec2 derivative(double t) const = 0;
    virtual Vec2 second_derivative(double t) const = 0;

    virtual double speed(double t) const;

    // T is defined by integral_t^(t+T) |Gamma'(u)| du = L.
    virtual double lookahead_increment(double t, double length, ArcLengthMethod method) const;
};

class VerticalLineCurve final : public Curve2D {
public:
    std::string name() const override;
    Vec2 position(double t) const override;
    Vec2 derivative(double t) const override;
    Vec2 second_derivative(double t) const override;
    double speed(double t) const override;
    double lookahead_increment(double t, double length, ArcLengthMethod method) const override;
};

class EllipseCurve final : public Curve2D {
public:
    EllipseCurve(double a_value, double b_value);

    std::string name() const override;
    Vec2 position(double t) const override;
    Vec2 derivative(double t) const override;
    Vec2 second_derivative(double t) const override;
    double speed(double t) const override;
    double lookahead_increment(double t, double length, ArcLengthMethod method) const override;

    double a() const { return a_; }
    double b() const { return b_; }

private:
    double numerical_lookahead(double t, double length) const;
    double arc_length_from(double t, double increment) const;

    double a_;
    double b_;
};

std::unique_ptr<Curve2D> make_curve(const std::string& curve_type, double ellipse_a, double ellipse_b);

}  // namespace gauss_trace
