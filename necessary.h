#pragma once

#include "curve.hpp"
#include "simulation.hpp"
#include "vec2.hpp"

// Backward-compatible alias for the original helper name.
inline double wrap(double x) {
    return gauss_trace::wrap_angle(x);
}
