#pragma once

// Ray + intersection queries. Ray direction need not be unit length, but
// distances returned by the *Hit functions are in units of `direction`.
// For Euclidean t, pass a unit direction.

#include "bromath/aabb.h"
#include "bromath/plane.h"
#include "bromath/scalar.h"
#include "bromath/sphere.h"
#include "bromath/vec.h"

#include <cmath>
#include <limits>

namespace bromath {

struct Ray {
    Vec3 origin    = {0, 0, 0};
    Vec3 direction = {0, 0, 1};
};

struct RayHit {
    bool  hit = false;
    float t   = 0.0f;
    Vec3  position{};
    Vec3  normal{};
};

inline constexpr Vec3 rat(const Ray& r, float t) {
    return r.origin + r.direction * t;
}

// Slab method (Williams). Returns t of first hit in [0, +inf), or hit=false.
// From outside, that is the entry point with the entered face's outward
// normal; from inside, the exit point with the exited face's outward normal.
// A zero direction hits nothing.
inline RayHit rIntersectAABB(const Ray& r, const AABB3& a) {
    float tmin = -std::numeric_limits<float>::infinity();
    float tmax =  std::numeric_limits<float>::infinity();
    const float* o = &r.origin.x;
    const float* d = &r.direction.x;
    const float* lo = &a.min.x;
    const float* hi = &a.max.x;
    int entryAxis = -1, exitAxis = -1;
    float entrySign = 1.0f, exitSign = 1.0f;
    for (int i = 0; i < 3; ++i) {
        if (std::fabs(d[i]) < 1e-20f) {
            if (o[i] < lo[i] || o[i] > hi[i]) return {};
            continue;
        }
        float inv = 1.0f / d[i];
        float t1 = (lo[i] - o[i]) * inv;
        float t2 = (hi[i] - o[i]) * inv;
        // Moving toward +axis the ray enters through the min face (outward
        // normal -1) and leaves through the max face (+1); toward -axis the
        // other way round.
        float lower = t1, upper = t2;
        float inSign = -1.0f, outSign = 1.0f;
        if (t1 > t2) { lower = t2; upper = t1; inSign = 1.0f; outSign = -1.0f; }
        if (lower > tmin) { tmin = lower; entryAxis = i; entrySign = inSign; }
        if (upper < tmax) { tmax = upper; exitAxis = i; exitSign = outSign; }
        if (tmin > tmax) return {};
    }
    if (exitAxis < 0) return {};  // zero direction: never reaches a face
    if (tmax < 0.0f) return {};
    const bool entering = tmin >= 0.0f;
    float t = entering ? tmin : tmax;
    Vec3 n{0, 0, 0};
    (&n.x)[entering ? entryAxis : exitAxis] = entering ? entrySign : exitSign;
    return { true, t, rat(r, t), n };
}

// Quadratic solution. Returns nearest non-negative t.
inline RayHit rIntersectSphere(const Ray& r, const Sphere& s) {
    Vec3 oc = r.origin - s.center;
    float a = vdot(r.direction, r.direction);
    if (!(a > 1e-30f)) return {};  // zero (or NaN) direction
    float b = 2.0f * vdot(oc, r.direction);
    float c = vdot(oc, oc) - s.radius * s.radius;
    float disc = b*b - 4*a*c;
    if (disc < 0.0f) return {};
    float sq = std::sqrt(disc);
    float t0 = (-b - sq) / (2*a);
    float t1 = (-b + sq) / (2*a);
    float t = t0 >= 0.0f ? t0 : t1;
    if (t < 0.0f) return {};
    Vec3 p = rat(r, t);
    return { true, t, p, vnorm(p - s.center) };
}

inline RayHit rIntersectPlane(const Ray& r, const Plane& p) {
    float denom = vdot(p.normal, r.direction);
    if (std::fabs(denom) < 1e-20f) return {};
    float t = -(vdot(p.normal, r.origin) + p.d) / denom;
    if (t < 0.0f) return {};
    return { true, t, rat(r, t), p.normal };
}

// Möller-Trumbore triangle intersection. Single-sided when backfaceCull = true.
inline RayHit rIntersectTriangle(const Ray& r, Vec3 v0, Vec3 v1, Vec3 v2,
                                 bool backfaceCull = false) {
    Vec3 e1 = v1 - v0;
    Vec3 e2 = v2 - v0;
    Vec3 pvec = vcross(r.direction, e2);
    float det = vdot(e1, pvec);
    if (backfaceCull) {
        if (det < 1e-20f) return {};
    } else {
        if (std::fabs(det) < 1e-20f) return {};
    }
    float invDet = 1.0f / det;
    Vec3 tvec = r.origin - v0;
    float u = vdot(tvec, pvec) * invDet;
    if (u < 0.0f || u > 1.0f) return {};
    Vec3 qvec = vcross(tvec, e1);
    float v = vdot(r.direction, qvec) * invDet;
    if (v < 0.0f || u + v > 1.0f) return {};
    float t = vdot(e2, qvec) * invDet;
    if (t < 0.0f) return {};
    Vec3 n = vnorm(vcross(e1, e2));
    return { true, t, rat(r, t), n };
}

} // namespace bromath
