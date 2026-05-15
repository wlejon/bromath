#pragma once

// Color: linear-float RGBA primary. Color8 is sRGB 8-bit-per-channel for
// pixel buffers / hex. Linear is the math-friendly form; conversion happens
// at the storage boundary.

#include "bromath/scalar.h"

#include <cmath>
#include <cstdint>

namespace bromath {

struct Color {
    float r = 0, g = 0, b = 0, a = 1;
    constexpr Color() = default;
    constexpr Color(float rr, float gg, float bb, float aa = 1.0f)
        : r(rr), g(gg), b(bb), a(aa) {}
};

struct Color8 {
    uint8_t r = 0, g = 0, b = 0, a = 255;
    constexpr Color8() = default;
    constexpr Color8(uint8_t rr, uint8_t gg, uint8_t bb, uint8_t aa = 255)
        : r(rr), g(gg), b(bb), a(aa) {}
};

inline constexpr Color clinear(float r, float g, float b, float a = 1.0f) {
    return {r, g, b, a};
}

inline constexpr Color clerp(Color a, Color b, float t) {
    return {
        lerp(a.r, b.r, t), lerp(a.g, b.g, t),
        lerp(a.b, b.b, t), lerp(a.a, b.a, t)
    };
}

// Gamma 2.2 approximation — fast, perceptually close to true sRGB.
inline float csrgbToLinear(float c) {
    if (c <= 0.04045f) return c / 12.92f;
    return std::pow((c + 0.055f) / 1.055f, 2.4f);
}

inline float clinearToSrgb(float c) {
    if (c <= 0.0031308f) return c * 12.92f;
    return 1.055f * std::pow(c, 1.0f / 2.4f) - 0.055f;
}

// Convert sRGB Color8 to linear Color (alpha is linear in both).
inline Color cfromColor8(Color8 c) {
    return {
        csrgbToLinear(c.r / 255.0f),
        csrgbToLinear(c.g / 255.0f),
        csrgbToLinear(c.b / 255.0f),
        c.a / 255.0f
    };
}

inline Color8 ctoColor8(Color c) {
    auto enc = [](float x) -> uint8_t {
        float s = clinearToSrgb(saturate(x)) * 255.0f + 0.5f;
        return (uint8_t)clamp(s, 0.0f, 255.0f);
    };
    return { enc(c.r), enc(c.g), enc(c.b),
             (uint8_t)clamp(c.a * 255.0f + 0.5f, 0.0f, 255.0f) };
}

// Parse "#RRGGBB" or "#RRGGBBAA" as sRGB and return linear Color. Returns
// fully transparent black on malformed input.
inline Color cfromHex(const char* hex) {
    if (!hex || hex[0] != '#') return {0, 0, 0, 0};
    auto hd = [](char c) -> int {
        if (c >= '0' && c <= '9') return c - '0';
        if (c >= 'a' && c <= 'f') return 10 + (c - 'a');
        if (c >= 'A' && c <= 'F') return 10 + (c - 'A');
        return -1;
    };
    int v[8] = {-1,-1,-1,-1,-1,-1,-1,-1};
    int n = 0;
    for (int i = 1; hex[i] && n < 8; ++i, ++n) v[n] = hd(hex[i]);
    if (n != 6 && n != 8) return {0, 0, 0, 0};
    for (int i = 0; i < n; ++i) if (v[i] < 0) return {0, 0, 0, 0};
    Color8 c{
        (uint8_t)(v[0]*16 + v[1]),
        (uint8_t)(v[2]*16 + v[3]),
        (uint8_t)(v[4]*16 + v[5]),
        (uint8_t)(n == 8 ? v[6]*16 + v[7] : 255)
    };
    return cfromColor8(c);
}

// HSV: h in [0,360) degrees, s and v in [0,1]. Returns linear Color (treats
// the HSV color as sRGB and converts).
inline Color cfromHSV(float h, float s, float v, float a = 1.0f) {
    h = std::fmod(h, 360.0f);
    if (h < 0.0f) h += 360.0f;
    float c = v * s;
    float x = c * (1.0f - std::fabs(std::fmod(h / 60.0f, 2.0f) - 1.0f));
    float m = v - c;
    float r=0, g=0, b=0;
    if      (h <  60.0f) { r=c; g=x; b=0; }
    else if (h < 120.0f) { r=x; g=c; b=0; }
    else if (h < 180.0f) { r=0; g=c; b=x; }
    else if (h < 240.0f) { r=0; g=x; b=c; }
    else if (h < 300.0f) { r=x; g=0; b=c; }
    else                 { r=c; g=0; b=x; }
    return {
        csrgbToLinear(r + m), csrgbToLinear(g + m), csrgbToLinear(b + m), a
    };
}

} // namespace bromath
