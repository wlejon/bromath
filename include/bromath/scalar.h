#pragma once

// Scalar utilities and constants. Header-only, POD-friendly, no dependencies
// beyond <cmath>. All functions are constexpr where the math permits.

#include <cmath>

namespace bromath {

inline constexpr float PI       = 3.14159265358979323846f;
inline constexpr float TWO_PI   = 6.28318530717958647692f;
inline constexpr float HALF_PI  = 1.57079632679489661923f;
inline constexpr float INV_PI   = 0.31830988618379067154f;
inline constexpr float INV_TWO_PI = 0.15915494309189533577f;
inline constexpr float DEG2RAD  = 0.01745329251994329577f;  // PI / 180
inline constexpr float RAD2DEG  = 57.2957795130823208768f;  // 180 / PI

inline constexpr float deg2rad(float d) { return d * DEG2RAD; }
inline constexpr float rad2deg(float r) { return r * RAD2DEG; }

template <typename T>
inline constexpr T min(T a, T b) { return a < b ? a : b; }
template <typename T>
inline constexpr T max(T a, T b) { return a > b ? a : b; }

template <typename T>
inline constexpr T clamp(T x, T lo, T hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

inline constexpr float saturate(float x) { return clamp(x, 0.0f, 1.0f); }

inline constexpr float lerp(float a, float b, float t) {
    return a + (b - a) * t;
}

// Alias matching GLSL/HLSL naming.
inline constexpr float mix(float a, float b, float t) {
    return a + (b - a) * t;
}

// Inverse lerp: where in [a,b] does x lie? Returns 0 when a==b.
inline constexpr float invLerp(float a, float b, float x) {
    float d = b - a;
    return d == 0.0f ? 0.0f : (x - a) / d;
}

// Remap x from [inMin,inMax] to [outMin,outMax].
inline constexpr float remap(float x, float inMin, float inMax,
                             float outMin, float outMax) {
    return lerp(outMin, outMax, invLerp(inMin, inMax, x));
}

// Hermite S-curve on [0,1]: 3t^2 - 2t^3, clamped at the edges.
inline constexpr float smoothstep(float edge0, float edge1, float x) {
    float t = saturate(invLerp(edge0, edge1, x));
    return t * t * (3.0f - 2.0f * t);
}

// 0-1 variant — same S-curve assuming x already in [0,1].
inline constexpr float smoothstep01(float x) {
    float t = saturate(x);
    return t * t * (3.0f - 2.0f * t);
}

// Quintic smootherstep: 6t^5 - 15t^4 + 10t^3 (Perlin's improved variant).
inline constexpr float smootherstep(float edge0, float edge1, float x) {
    float t = saturate(invLerp(edge0, edge1, x));
    return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
}

inline constexpr float step(float edge, float x) {
    return x < edge ? 0.0f : 1.0f;
}

template <typename T>
inline constexpr T sign(T x) {
    return T((x > T(0)) - (x < T(0)));
}

template <typename T>
inline constexpr T abs(T x) { return x < T(0) ? -x : x; }

// Square — handy enough to be worth a name.
template <typename T>
inline constexpr T sqr(T x) { return x * x; }

// Fast approximate equality for floats.
inline bool nearlyEqual(float a, float b, float eps = 1e-5f) {
    return std::fabs(a - b) <= eps;
}

} // namespace bromath
