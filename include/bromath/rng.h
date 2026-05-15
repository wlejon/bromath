#pragma once

// Deterministic RNG. SplitMix64 is the engine — caller threads a uint64_t
// state through, no global mutation. Same primitive previously vendored
// inside broflora; consolidated here so all siblings reach for the same
// generator and get reproducible results from a shared seed.

#include "bromath/scalar.h"
#include "bromath/vec.h"

#include <cmath>
#include <cstdint>

namespace bromath {

inline uint64_t splitmix64(uint64_t& state) {
    uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// Uniform float in [0, 1). Uses the top 24 bits as a single-precision mantissa.
inline float randFloat01(uint64_t& state) {
    return (splitmix64(state) >> 40) * (1.0f / 16777216.0f);
}

inline float randSigned(uint64_t& state) {
    return randFloat01(state) * 2.0f - 1.0f;
}

inline float randRange(uint64_t& state, float lo, float hi) {
    return lerp(lo, hi, randFloat01(state));
}

inline int randInt(uint64_t& state, int loInclusive, int hiInclusive) {
    if (hiInclusive < loInclusive) return loInclusive;
    uint64_t span = (uint64_t)(hiInclusive - loInclusive) + 1ULL;
    return loInclusive + (int)(splitmix64(state) % span);
}

// Box-Muller standard normal sample.
inline float randNormal(uint64_t& state) {
    float u1 = randFloat01(state);
    if (u1 < 1e-7f) u1 = 1e-7f;
    float u2 = randFloat01(state);
    return std::sqrt(-2.0f * std::log(u1)) * std::cos(TWO_PI * u2);
}

inline Vec2 randGaussian2D(uint64_t& state, float sigma) {
    return { randNormal(state) * sigma, randNormal(state) * sigma };
}

// Uniformly distributed point inside the unit disc (XY plane).
inline Vec2 randInUnitDisc(uint64_t& state) {
    // Rejection sampling — simple and uniform.
    for (int i = 0; i < 32; ++i) {
        Vec2 p{ randSigned(state), randSigned(state) };
        if (vlen2(p) <= 1.0f) return p;
    }
    return {0, 0};
}

// Uniformly distributed point inside the unit sphere.
inline Vec3 randInUnitSphere(uint64_t& state) {
    for (int i = 0; i < 32; ++i) {
        Vec3 p{ randSigned(state), randSigned(state), randSigned(state) };
        if (vlen2(p) <= 1.0f) return p;
    }
    return {0, 0, 0};
}

// Uniformly distributed point on the unit sphere surface.
inline Vec3 randOnUnitSphere(uint64_t& state) {
    float z = randSigned(state);
    float a = randFloat01(state) * TWO_PI;
    float r = std::sqrt(1.0f - z * z);
    return { r * std::cos(a), r * std::sin(a), z };
}

} // namespace bromath
