#pragma once

// One-pole exponential smoother for parameter ramps. Lifted from broaudio's
// per-voice gain/pan smoothing — useful any time you want to chase a target
// without zipper artifacts: scene parameter automation, camera follow,
// physics-decoupled UI animation.
//
// Coefficient is chosen so the smoothed value reaches ~95% of the target in
// `timeMs` milliseconds at the given sample rate (or tick rate).

#include "bromath/scalar.h"

#include <cmath>

namespace bromath {

struct Smoother {
    float current = 0.0f;
    float target  = 0.0f;
    float coeff   = 0.0f;
};

// Set the smoothing time constant. After `timeMs` of ticking, the smoother
// will have closed ~95% of the gap between current and target.
inline void smootherSetTime(Smoother& s, float timeMs, float sampleRate) {
    if (timeMs <= 0.0f || sampleRate <= 0.0f) { s.coeff = 1.0f; return; }
    float samples = timeMs * 0.001f * sampleRate;
    s.coeff = 1.0f - std::pow(0.05f, 1.0f / samples);
}

inline void smootherReset(Smoother& s, float value) {
    s.current = s.target = value;
}

inline void smootherTarget(Smoother& s, float t) { s.target = t; }

inline float smootherTick(Smoother& s) {
    s.current += (s.target - s.current) * s.coeff;
    return s.current;
}

// Advance N ticks at once — useful when the host loop is block-based.
// Closed form: each tick keeps (1 - coeff) of the gap, so N ticks keep
// (1 - coeff)^N of it, and the cost does not grow with N.
inline float smootherTickN(Smoother& s, int n) {
    if (n <= 0) return s.current;
    double keep = std::pow(1.0 - static_cast<double>(s.coeff), static_cast<double>(n));
    s.current = static_cast<float>(s.target + (static_cast<double>(s.current) - s.target) * keep);
    return s.current;
}

} // namespace bromath
