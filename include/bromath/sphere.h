#pragma once

// Sphere primitive. center + radius; queries for containment, intersection,
// and a closed-form sphere-sphere intersection volume (lifted from broflora).

#include "bromath/scalar.h"
#include "bromath/vec.h"

#include <cmath>

namespace bromath {

struct Sphere {
    Vec3 center = {0, 0, 0};
    float radius = 0.0f;
};

inline constexpr bool scontains(const Sphere& s, Vec3 p) {
    return vdist2(s.center, p) <= s.radius * s.radius;
}

inline constexpr bool sintersects(const Sphere& a, const Sphere& b) {
    float r = a.radius + b.radius;
    return vdist2(a.center, b.center) <= r * r;
}

// Volume of intersection of two spheres. Returns 0 when disjoint; returns
// the smaller sphere's volume when one contains the other.
inline float sintersectVolume(const Sphere& a, const Sphere& b) {
    float d = vdist(a.center, b.center);
    if (d >= a.radius + b.radius) return 0.0f;
    if (d + min(a.radius, b.radius) <= max(a.radius, b.radius)) {
        float r = min(a.radius, b.radius);
        return (4.0f / 3.0f) * PI * r * r * r;
    }
    // Lens-volume formula: sum of two spherical caps.
    float r1 = a.radius, r2 = b.radius;
    float t = (PI * sqr(r1 + r2 - d)
               * (d*d + 2*d*r1 - 3*r1*r1 + 2*d*r2 + 6*r1*r2 - 3*r2*r2))
              / (12.0f * d);
    return t;
}

} // namespace bromath
