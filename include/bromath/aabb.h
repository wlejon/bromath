#pragma once

// Axis-aligned bounding boxes — 2D (AABB2) and 3D (AABB3). Stored as
// min/max corners. Empty boxes have min > max along some axis; build
// new empties with aempty() and merge points/boxes to grow.

#include "bromath/mat.h"
#include "bromath/scalar.h"
#include "bromath/vec.h"

#include <limits>

namespace bromath {

struct AABB2 {
    Vec2 min = {0, 0};
    Vec2 max = {0, 0};
};

struct AABB3 {
    Vec3 min = {0, 0, 0};
    Vec3 max = {0, 0, 0};
};

inline AABB3 aempty3() {
    float inf = std::numeric_limits<float>::infinity();
    return { { inf, inf, inf }, { -inf, -inf, -inf } };
}

inline AABB2 aempty2() {
    float inf = std::numeric_limits<float>::infinity();
    return { { inf, inf }, { -inf, -inf } };
}

inline constexpr bool aisEmpty(const AABB3& a) {
    return a.min.x > a.max.x || a.min.y > a.max.y || a.min.z > a.max.z;
}

inline constexpr Vec3 acenter(const AABB3& a) {
    return (a.min + a.max) * 0.5f;
}

inline constexpr Vec3 aextent(const AABB3& a) { return a.max - a.min; }
inline constexpr Vec3 ahalfExtent(const AABB3& a) { return aextent(a) * 0.5f; }

inline constexpr bool acontains(const AABB3& a, Vec3 p) {
    return p.x >= a.min.x && p.x <= a.max.x
        && p.y >= a.min.y && p.y <= a.max.y
        && p.z >= a.min.z && p.z <= a.max.z;
}

inline constexpr bool aintersects(const AABB3& a, const AABB3& b) {
    return a.min.x <= b.max.x && a.max.x >= b.min.x
        && a.min.y <= b.max.y && a.max.y >= b.min.y
        && a.min.z <= b.max.z && a.max.z >= b.min.z;
}

inline constexpr AABB3 aexpand(const AABB3& a, Vec3 p) {
    return { vmin(a.min, p), vmax(a.max, p) };
}

inline constexpr AABB3 amerge(const AABB3& a, const AABB3& b) {
    return { vmin(a.min, b.min), vmax(a.max, b.max) };
}

// Smallest AABB containing all `count` points (stride in floats, default 3).
inline AABB3 afromPoints(const float* pts, int count, int stride = 3) {
    if (count <= 0) return aempty3();
    AABB3 a;
    a.min = a.max = { pts[0], pts[1], pts[2] };
    for (int i = 1; i < count; ++i) {
        const float* p = pts + i * stride;
        a = aexpand(a, Vec3{p[0], p[1], p[2]});
    }
    return a;
}

// Transform an AABB by a 4x4 affine matrix and return the axis-aligned
// bound of the result. Uses Arvo's trick: O(1) per corner via row sums.
inline AABB3 atransform(const AABB3& a, const Mat4& m) {
    Vec3 c = acenter(a);
    Vec3 e = ahalfExtent(a);
    Vec3 newC = mtransformPoint(m, c);
    Vec3 newE{
        std::fabs(m.at(0,0)) * e.x + std::fabs(m.at(0,1)) * e.y + std::fabs(m.at(0,2)) * e.z,
        std::fabs(m.at(1,0)) * e.x + std::fabs(m.at(1,1)) * e.y + std::fabs(m.at(1,2)) * e.z,
        std::fabs(m.at(2,0)) * e.x + std::fabs(m.at(2,1)) * e.y + std::fabs(m.at(2,2)) * e.z
    };
    return { newC - newE, newC + newE };
}

// --- AABB2 ---
inline constexpr bool acontains(const AABB2& a, Vec2 p) {
    return p.x >= a.min.x && p.x <= a.max.x
        && p.y >= a.min.y && p.y <= a.max.y;
}
inline constexpr bool aintersects(const AABB2& a, const AABB2& b) {
    return a.min.x <= b.max.x && a.max.x >= b.min.x
        && a.min.y <= b.max.y && a.max.y >= b.min.y;
}

} // namespace bromath
