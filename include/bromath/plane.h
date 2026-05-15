#pragma once

// Plane in implicit form: dot(normal, p) + d = 0. Normal is expected to
// be unit-length for distance queries to return Euclidean values.

#include "bromath/vec.h"

namespace bromath {

struct Plane {
    Vec3 normal = {0, 1, 0};
    float d = 0.0f;
};

inline Plane pfromPointNormal(Vec3 point, Vec3 normal) {
    Vec3 n = vnorm(normal);
    return { n, -vdot(n, point) };
}

inline Plane pfromPoints(Vec3 a, Vec3 b, Vec3 c) {
    Vec3 n = vnorm(vcross(b - a, c - a));
    return { n, -vdot(n, a) };
}

// Signed distance from point to plane. Positive on the side `normal` points to.
inline constexpr float psignedDistance(const Plane& p, Vec3 q) {
    return vdot(p.normal, q) + p.d;
}

// Project point `q` onto the plane.
inline constexpr Vec3 pproject(const Plane& p, Vec3 q) {
    return q - p.normal * psignedDistance(p, q);
}

} // namespace bromath
