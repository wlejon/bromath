#pragma once

// Capsule primitive + segment/segment closest distance.
//
// A capsule is a line segment (p0..p1) inflated by `radius` — the natural
// bounding volume for a branch, limb, or any swept sphere. Unlike a sphere,
// a capsule is *directional*: two capsules whose cores cross interpenetrate
// even when their midpoints are far apart, which a sphere-vs-sphere test
// misses. Used by collision code that needs to reason about elongated
// geometry (e.g. branches growing through one another).

#include "bromath/scalar.h"
#include "bromath/vec.h"

#include <cmath>

namespace bromath {

struct Capsule {
    Vec3  p0     = {0, 0, 0};
    Vec3  p1     = {0, 0, 0};
    float radius = 0.0f;
};

// Squared closest distance between segments [p1,q1] and [p2,q2], with the
// clamped parameters s,t of the closest points (c1 = p1 + s*(q1-p1),
// c2 = p2 + t*(q2-p2)). Ericson, Real-Time Collision Detection §5.1.9.
// Handles degenerate (zero-length) segments by collapsing to point tests.
inline float closestSegmentSegment2(Vec3 p1, Vec3 q1, Vec3 p2, Vec3 q2,
                                    float& s, float& t) {
    const float EPS = 1e-12f;
    Vec3 d1 = q1 - p1;   // direction of segment S1
    Vec3 d2 = q2 - p2;   // direction of segment S2
    Vec3 r  = p1 - p2;
    float a = vdot(d1, d1);   // squared length of S1
    float e = vdot(d2, d2);   // squared length of S2
    float f = vdot(d2, r);

    if (a <= EPS && e <= EPS) {
        // Both segments degenerate to points.
        s = 0.0f; t = 0.0f;
        Vec3 c = p1 - p2;
        return vdot(c, c);
    }
    if (a <= EPS) {
        // First segment degenerate.
        s = 0.0f;
        t = clamp(f / e, 0.0f, 1.0f);
    } else {
        float c = vdot(d1, r);
        if (e <= EPS) {
            // Second segment degenerate.
            t = 0.0f;
            s = clamp(-c / a, 0.0f, 1.0f);
        } else {
            // General non-degenerate case.
            float b = vdot(d1, d2);
            float denom = a * e - b * b;   // always >= 0
            s = (denom > EPS) ? clamp((b * f - c * e) / denom, 0.0f, 1.0f) : 0.0f;
            t = (b * s + f) / e;
            // If t out of [0,1], clamp and recompute s for the new t.
            if (t < 0.0f) {
                t = 0.0f;
                s = clamp(-c / a, 0.0f, 1.0f);
            } else if (t > 1.0f) {
                t = 1.0f;
                s = clamp((b - c) / a, 0.0f, 1.0f);
            }
        }
    }
    Vec3 c1 = p1 + d1 * s;
    Vec3 c2 = p2 + d2 * t;
    Vec3 cc = c1 - c2;
    return vdot(cc, cc);
}

inline float segmentSegmentDistance(Vec3 p1, Vec3 q1, Vec3 p2, Vec3 q2) {
    float s, t;
    return std::sqrt(closestSegmentSegment2(p1, q1, p2, q2, s, t));
}

// Interpenetration depth of two capsules: how far their surfaces overlap,
// 0 when disjoint. = max(0, r0 + r1 - distance between their core segments).
inline float capsulePenetration(const Capsule& a, const Capsule& b) {
    float s, t;
    float d2 = closestSegmentSegment2(a.p0, a.p1, b.p0, b.p1, s, t);
    float rr = a.radius + b.radius;
    if (d2 >= rr * rr) return 0.0f;
    return rr - std::sqrt(d2);
}

inline bool capsulesIntersect(const Capsule& a, const Capsule& b) {
    float s, t;
    float d2 = closestSegmentSegment2(a.p0, a.p1, b.p0, b.p1, s, t);
    float rr = a.radius + b.radius;
    return d2 <= rr * rr;
}

} // namespace bromath
