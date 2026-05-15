#pragma once

// Vec2 / Vec3 with free-function ops. POD aggregates — trivially copyable,
// bindings-friendly. Free functions (vdot, vcross, vnorm, ...) compose
// naturally with raw float arrays the way bromesh stores mesh data.

#include "bromath/scalar.h"

#include <cmath>

namespace bromath {

struct Vec2 {
    float x = 0, y = 0;
    constexpr Vec2() = default;
    constexpr Vec2(float xx, float yy) : x(xx), y(yy) {}
};

struct Vec3 {
    float x = 0, y = 0, z = 0;
    constexpr Vec3() = default;
    constexpr Vec3(float xx, float yy, float zz) : x(xx), y(yy), z(zz) {}
};

// --- Vec2 ---
inline constexpr Vec2 operator+(Vec2 a, Vec2 b) { return {a.x+b.x, a.y+b.y}; }
inline constexpr Vec2 operator-(Vec2 a, Vec2 b) { return {a.x-b.x, a.y-b.y}; }
inline constexpr Vec2 operator-(Vec2 a)         { return {-a.x, -a.y}; }
inline constexpr Vec2 operator*(Vec2 a, float s){ return {a.x*s, a.y*s}; }
inline constexpr Vec2 operator*(float s, Vec2 a){ return {a.x*s, a.y*s}; }
inline constexpr Vec2 operator/(Vec2 a, float s){ return {a.x/s, a.y/s}; }
inline constexpr Vec2& operator+=(Vec2& a, Vec2 b){ a.x+=b.x; a.y+=b.y; return a; }
inline constexpr Vec2& operator-=(Vec2& a, Vec2 b){ a.x-=b.x; a.y-=b.y; return a; }
inline constexpr Vec2& operator*=(Vec2& a, float s){ a.x*=s; a.y*=s; return a; }

inline constexpr float vdot(Vec2 a, Vec2 b) { return a.x*b.x + a.y*b.y; }
inline constexpr float vlen2(Vec2 a)        { return vdot(a, a); }
inline float           vlen(Vec2 a)         { return std::sqrt(vlen2(a)); }
inline float           vdist(Vec2 a, Vec2 b){ return vlen(a - b); }
inline constexpr float vdist2(Vec2 a, Vec2 b){ Vec2 d = a - b; return vdot(d, d); }
// 2D scalar cross — z-component of (a,0) x (b,0). Useful for winding/turn tests.
inline constexpr float vcross(Vec2 a, Vec2 b){ return a.x*b.y - a.y*b.x; }

inline Vec2 vnorm(Vec2 a) {
    float L = vlen(a);
    return (L > 1e-20f) ? Vec2{a.x/L, a.y/L} : Vec2{0, 0};
}
inline Vec2 vnormOr(Vec2 a, Vec2 fallback) {
    float L = vlen(a);
    return (L > 1e-20f) ? Vec2{a.x/L, a.y/L} : fallback;
}
inline constexpr Vec2 vlerp(Vec2 a, Vec2 b, float t) {
    return {lerp(a.x, b.x, t), lerp(a.y, b.y, t)};
}

// --- Vec3 ---
inline constexpr Vec3 operator+(Vec3 a, Vec3 b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
inline constexpr Vec3 operator-(Vec3 a, Vec3 b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
inline constexpr Vec3 operator-(Vec3 a)         { return {-a.x, -a.y, -a.z}; }
inline constexpr Vec3 operator*(Vec3 a, float s){ return {a.x*s, a.y*s, a.z*s}; }
inline constexpr Vec3 operator*(float s, Vec3 a){ return {a.x*s, a.y*s, a.z*s}; }
inline constexpr Vec3 operator/(Vec3 a, float s){ return {a.x/s, a.y/s, a.z/s}; }
inline constexpr Vec3& operator+=(Vec3& a, Vec3 b){ a.x+=b.x; a.y+=b.y; a.z+=b.z; return a; }
inline constexpr Vec3& operator-=(Vec3& a, Vec3 b){ a.x-=b.x; a.y-=b.y; a.z-=b.z; return a; }
inline constexpr Vec3& operator*=(Vec3& a, float s){ a.x*=s; a.y*=s; a.z*=s; return a; }

inline constexpr float vdot(Vec3 a, Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline constexpr Vec3  vcross(Vec3 a, Vec3 b) {
    return { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x };
}
inline constexpr float vlen2(Vec3 a)        { return vdot(a, a); }
inline float           vlen(Vec3 a)         { return std::sqrt(vlen2(a)); }
inline float           vdist(Vec3 a, Vec3 b){ return vlen(a - b); }
inline constexpr float vdist2(Vec3 a, Vec3 b){ Vec3 d = a - b; return vdot(d, d); }

inline Vec3 vnorm(Vec3 a) {
    float L = vlen(a);
    return (L > 1e-20f) ? Vec3{a.x/L, a.y/L, a.z/L} : Vec3{0, 0, 0};
}
inline Vec3 vnormOr(Vec3 a, Vec3 fallback) {
    float L = vlen(a);
    return (L > 1e-20f) ? Vec3{a.x/L, a.y/L, a.z/L} : fallback;
}
inline constexpr Vec3 vlerp(Vec3 a, Vec3 b, float t) {
    return {lerp(a.x, b.x, t), lerp(a.y, b.y, t), lerp(a.z, b.z, t)};
}

// Reflect `i` about unit normal `n`. n must be unit length.
inline constexpr Vec3 vreflect(Vec3 i, Vec3 n) {
    return i - n * (2.0f * vdot(i, n));
}

// Project `a` onto `b`. Returns {0,0,0} if b is degenerate.
inline Vec3 vproject(Vec3 a, Vec3 b) {
    float d = vdot(b, b);
    if (d < 1e-20f) return {0, 0, 0};
    return b * (vdot(a, b) / d);
}

// Component-wise min/max — useful for AABB construction.
inline constexpr Vec3 vmin(Vec3 a, Vec3 b) {
    return { min(a.x,b.x), min(a.y,b.y), min(a.z,b.z) };
}
inline constexpr Vec3 vmax(Vec3 a, Vec3 b) {
    return { max(a.x,b.x), max(a.y,b.y), max(a.z,b.z) };
}

// Pick an arbitrary unit vector perpendicular to `n`. n need not be unit.
inline Vec3 vperpendicular(Vec3 n) {
    Vec3 axis = (std::fabs(n.x) < 0.9f) ? Vec3{1,0,0} : Vec3{0,1,0};
    return vnorm(vcross(n, axis));
}

} // namespace bromath
