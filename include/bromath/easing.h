#pragma once

// Easing functions (Penner set). All map progress t in [0,1] to eased
// progress with f(0) = 0 and f(1) = 1. back/elastic intentionally overshoot
// outside [0,1] mid-curve; bounce stays inside. Inputs outside [0,1] are the
// caller's responsibility (clamp first for hard stops, or let back/elastic
// extrapolate for anticipation/follow-through).

#include "bromath/scalar.h"

#include <cmath>

namespace bromath {

inline float easeLinear(float t) { return t; }

// --- Power curves -----------------------------------------------------------

inline float easeQuadIn(float t)  { return t * t; }
inline float easeQuadOut(float t) { float u = 1.0f - t; return 1.0f - u * u; }
inline float easeQuadInOut(float t) {
    if (t < 0.5f) return 2.0f * t * t;
    float u = 1.0f - t; return 1.0f - 2.0f * u * u;
}

inline float easeCubicIn(float t)  { return t * t * t; }
inline float easeCubicOut(float t) { float u = 1.0f - t; return 1.0f - u * u * u; }
inline float easeCubicInOut(float t) {
    if (t < 0.5f) return 4.0f * t * t * t;
    float u = 1.0f - t; return 1.0f - 4.0f * u * u * u;
}

inline float easeQuartIn(float t)  { float s = t * t; return s * s; }
inline float easeQuartOut(float t) { float u = 1.0f - t; float s = u * u; return 1.0f - s * s; }
inline float easeQuartInOut(float t) {
    if (t < 0.5f) { float s = t * t; return 8.0f * s * s; }
    float u = 1.0f - t; float s = u * u; return 1.0f - 8.0f * s * s;
}

inline float easeQuintIn(float t)  { float s = t * t; return s * s * t; }
inline float easeQuintOut(float t) { float u = 1.0f - t; float s = u * u; return 1.0f - s * s * u; }
inline float easeQuintInOut(float t) {
    if (t < 0.5f) { float s = t * t; return 16.0f * s * s * t; }
    float u = 1.0f - t; float s = u * u; return 1.0f - 16.0f * s * s * u;
}

// --- Trig / exponential -----------------------------------------------------

inline float easeSineIn(float t)  { return 1.0f - std::cos(t * HALF_PI); }
inline float easeSineOut(float t) { return std::sin(t * HALF_PI); }
inline float easeSineInOut(float t) { return -(std::cos(PI * t) - 1.0f) * 0.5f; }

inline float easeExpoIn(float t) {
    return t <= 0.0f ? 0.0f : std::pow(2.0f, 10.0f * t - 10.0f);
}
inline float easeExpoOut(float t) {
    return t >= 1.0f ? 1.0f : 1.0f - std::pow(2.0f, -10.0f * t);
}
inline float easeExpoInOut(float t) {
    if (t <= 0.0f) return 0.0f;
    if (t >= 1.0f) return 1.0f;
    return t < 0.5f ? std::pow(2.0f, 20.0f * t - 10.0f) * 0.5f
                    : (2.0f - std::pow(2.0f, -20.0f * t + 10.0f)) * 0.5f;
}

inline float easeCircIn(float t)  { return 1.0f - std::sqrt(max(0.0f, 1.0f - t * t)); }
inline float easeCircOut(float t) { float u = t - 1.0f; return std::sqrt(max(0.0f, 1.0f - u * u)); }
inline float easeCircInOut(float t) {
    if (t < 0.5f) {
        float x = 2.0f * t;
        return (1.0f - std::sqrt(max(0.0f, 1.0f - x * x))) * 0.5f;
    }
    float x = -2.0f * t + 2.0f;
    return (std::sqrt(max(0.0f, 1.0f - x * x)) + 1.0f) * 0.5f;
}

// --- Overshoot / oscillation ------------------------------------------------

inline float easeBackIn(float t) {
    constexpr float c1 = 1.70158f, c3 = c1 + 1.0f;
    return c3 * t * t * t - c1 * t * t;
}
inline float easeBackOut(float t) {
    constexpr float c1 = 1.70158f, c3 = c1 + 1.0f;
    float u = t - 1.0f;
    return 1.0f + c3 * u * u * u + c1 * u * u;
}
inline float easeBackInOut(float t) {
    constexpr float c1 = 1.70158f, c2 = c1 * 1.525f;
    if (t < 0.5f) {
        float x = 2.0f * t;
        return (x * x * ((c2 + 1.0f) * x - c2)) * 0.5f;
    }
    float x = 2.0f * t - 2.0f;
    return (x * x * ((c2 + 1.0f) * x + c2) + 2.0f) * 0.5f;
}

inline float easeElasticIn(float t) {
    if (t <= 0.0f) return 0.0f;
    if (t >= 1.0f) return 1.0f;
    constexpr float c4 = TWO_PI / 3.0f;
    return -std::pow(2.0f, 10.0f * t - 10.0f) * std::sin((10.0f * t - 10.75f) * c4);
}
inline float easeElasticOut(float t) {
    if (t <= 0.0f) return 0.0f;
    if (t >= 1.0f) return 1.0f;
    constexpr float c4 = TWO_PI / 3.0f;
    return std::pow(2.0f, -10.0f * t) * std::sin((10.0f * t - 0.75f) * c4) + 1.0f;
}
inline float easeElasticInOut(float t) {
    if (t <= 0.0f) return 0.0f;
    if (t >= 1.0f) return 1.0f;
    constexpr float c5 = TWO_PI / 4.5f;
    return t < 0.5f
        ? -(std::pow(2.0f, 20.0f * t - 10.0f) * std::sin((20.0f * t - 11.125f) * c5)) * 0.5f
        : std::pow(2.0f, -20.0f * t + 10.0f) * std::sin((20.0f * t - 11.125f) * c5) * 0.5f + 1.0f;
}

inline float easeBounceOut(float t) {
    constexpr float n1 = 7.5625f, d1 = 2.75f;
    if (t < 1.0f / d1)  return n1 * t * t;
    if (t < 2.0f / d1)  { t -= 1.5f / d1;  return n1 * t * t + 0.75f; }
    if (t < 2.5f / d1)  { t -= 2.25f / d1; return n1 * t * t + 0.9375f; }
    t -= 2.625f / d1;   return n1 * t * t + 0.984375f;
}
inline float easeBounceIn(float t) { return 1.0f - easeBounceOut(1.0f - t); }
inline float easeBounceInOut(float t) {
    return t < 0.5f ? (1.0f - easeBounceOut(1.0f - 2.0f * t)) * 0.5f
                    : (1.0f + easeBounceOut(2.0f * t - 1.0f)) * 0.5f;
}

} // namespace bromath
