#pragma once

// Angle utilities. All angles in radians unless noted. Wrapping uses
// floor-based reduction to avoid an fmod libcall.

#include "bromath/scalar.h"

#include <cmath>

namespace bromath {

// Fold angle into [-PI, PI].
inline float wrapAngle(float a) {
    const float x = a + PI;
    return x - TWO_PI * std::floor(x * INV_TWO_PI) - PI;
}

// Fold angle into [0, TWO_PI).
inline float wrapAngle2Pi(float a) {
    return a - TWO_PI * std::floor(a * INV_TWO_PI);
}

// Shortest signed delta from `from` to `to`, in [-PI, PI].
inline float angleDelta(float from, float to) {
    return wrapAngle(to - from);
}

// Shortest-arc lerp between two angles (interpolates through the short way).
inline float angleLerp(float from, float to, float t) {
    return wrapAngle(from + angleDelta(from, to) * t);
}

} // namespace bromath
