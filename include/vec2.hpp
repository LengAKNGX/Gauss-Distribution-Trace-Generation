#pragma once

#include <cmath>
#include <stdexcept>

namespace gauss_trace {

struct Vec2 {
    double x = 0.0;
    double y = 0.0;

    Vec2() = default;
    Vec2(double x_value, double y_value) : x(x_value), y(y_value) {}

    Vec2 operator+(const Vec2& other) const { return {x + other.x, y + other.y}; }
    Vec2 operator-(const Vec2& other) const { return {x - other.x, y - other.y}; }
    Vec2 operator*(double scalar) const { return {x * scalar, y * scalar}; }
    Vec2 operator/(double scalar) const { return {x / scalar, y / scalar}; }

    Vec2& operator+=(const Vec2& other) {
        x += other.x;
        y += other.y;
        return *this;
    }
};

inline Vec2 operator*(double scalar, const Vec2& value) {
    return value * scalar;
}

inline double dot(const Vec2& lhs, const Vec2& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y;
}

inline double norm_squared(const Vec2& value) {
    return dot(value, value);
}

inline double norm(const Vec2& value) {
    return std::sqrt(norm_squared(value));
}

inline Vec2 normalized(const Vec2& value) {
    const double length = norm(value);
    if (length == 0.0) {
        throw std::runtime_error("cannot normalize zero vector");
    }
    return value / length;
}

}  // namespace gauss_trace
