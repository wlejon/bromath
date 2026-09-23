#pragma once

// Curve evaluators. All evaluate a single parameter t (typically in [0,1])
// to a Vec3. Float-valued variants used for CSS-style easing live here too.

#include "bromath/scalar.h"
#include "bromath/vec.h"

#include <cmath>

namespace bromath {

// Cubic Bézier in 1D — useful for easing curves (CSS cubic-bezier).
// Control points (0,p1y) and (1,p2y) are implicit at (p1x, p1y) / (p2x, p2y).
struct CubicEase {
    float p1x = 0.25f, p1y = 0.1f;
    float p2x = 0.25f, p2y = 1.0f;
};

// Sample y for given x in [0,1]. Solves t from x via Newton iteration, then
// evaluates the y polynomial.
inline float ccubicEase(const CubicEase& c, float x) {
    if (x <= 0.0f) return 0.0f;
    if (x >= 1.0f) return 1.0f;
    auto sampleX = [&](float t) {
        float u = 1.0f - t;
        return 3.0f*u*u*t*c.p1x + 3.0f*u*t*t*c.p2x + t*t*t;
    };
    auto sampleDX = [&](float t) {
        float u = 1.0f - t;
        return 3.0f*u*u*c.p1x + 6.0f*u*t*(c.p2x - c.p1x) + 3.0f*t*t*(1.0f - c.p2x);
    };
    float t = x;
    bool solved = false;
    for (int i = 0; i < 8; ++i) {
        float dx = sampleX(t) - x;
        if (std::fabs(dx) < 1e-6f) { solved = true; break; }
        float d = sampleDX(t);
        if (std::fabs(d) < 1e-6f) break;
        t -= dx / d;
        t = clamp(t, 0.0f, 1.0f);
    }
    if (!solved) {
        // Newton stalls where x(t) is flat (a control x of 0 or 1) or
        // oscillates on steep curves; x(t) is monotone for control x in
        // [0, 1], so bisect instead.
        float lo = 0.0f, hi = 1.0f;
        t = x;
        for (int i = 0; i < 40; ++i) {
            float sx = sampleX(t);
            if (std::fabs(sx - x) < 1e-7f) break;
            if (sx < x) lo = t; else hi = t;
            t = 0.5f * (lo + hi);
        }
    }
    float u = 1.0f - t;
    return 3.0f*u*u*t*c.p1y + 3.0f*u*t*t*c.p2y + t*t*t;
}

// 3D cubic Bézier with four explicit control points.
inline Vec3 cbezier(Vec3 p0, Vec3 p1, Vec3 p2, Vec3 p3, float t) {
    float u = 1.0f - t;
    float b0 = u*u*u;
    float b1 = 3.0f*u*u*t;
    float b2 = 3.0f*u*t*t;
    float b3 = t*t*t;
    return p0*b0 + p1*b1 + p2*b2 + p3*b3;
}

// Derivative of a cubic Bézier — useful for tangents.
inline Vec3 cbezierTangent(Vec3 p0, Vec3 p1, Vec3 p2, Vec3 p3, float t) {
    float u = 1.0f - t;
    return 3.0f*u*u*(p1 - p0) + 6.0f*u*t*(p2 - p1) + 3.0f*t*t*(p3 - p2);
}

// Cubic Hermite — used by glTF CubicSpline animation channels. m0 is the
// tangent leaving p0 and m1 the tangent arriving at p1 (for glTF: p0's
// out-tangent and p1's in-tangent, each scaled by the keyframe interval).
inline Vec3 chermite(Vec3 p0, Vec3 m0, Vec3 p1, Vec3 m1, float t) {
    float t2 = t * t;
    float t3 = t2 * t;
    float h00 =  2.0f*t3 - 3.0f*t2 + 1.0f;
    float h10 =       t3 - 2.0f*t2 + t;
    float h01 = -2.0f*t3 + 3.0f*t2;
    float h11 =       t3 -      t2;
    return p0*h00 + m0*h10 + p1*h01 + m1*h11;
}

// Centripetal Catmull-Rom (Lee parameter alpha=0.5). Passes through p1 and
// p2; p0 and p3 are the surrounding control points. Smoother than uniform
// Catmull-Rom when points are non-uniformly spaced.
inline Vec3 ccatmullRom(Vec3 p0, Vec3 p1, Vec3 p2, Vec3 p3, float t) {
    auto tj = [](float ti, Vec3 a, Vec3 b) {
        return std::pow(vdist(a, b), 0.5f) + ti;
    };
    // A segment of zero length has nothing to interpolate.
    if (vdist2(p1, p2) <= 1e-24f) return p1;
    // An end control point repeated onto its neighbour (the usual way to
    // start or end a spline) is replaced by the reflection of the other
    // neighbour, keeping the curve cubic instead of falling back to a lerp.
    if (vdist2(p0, p1) <= 1e-24f) p0 = p1 + (p1 - p2);
    if (vdist2(p3, p2) <= 1e-24f) p3 = p2 + (p2 - p1);
    float t0 = 0.0f;
    float t1 = tj(t0, p0, p1);
    float t2 = tj(t1, p1, p2);
    float t3 = tj(t2, p2, p3);
    if (!(t1 > t0 && t2 > t1 && t3 > t2)) {
        // Knots still coincide (non-finite input): plain lerp p1->p2.
        return vlerp(p1, p2, t);
    }
    float tt = lerp(t1, t2, t);
    Vec3 a1 = p0 * ((t1 - tt) / (t1 - t0)) + p1 * ((tt - t0) / (t1 - t0));
    Vec3 a2 = p1 * ((t2 - tt) / (t2 - t1)) + p2 * ((tt - t1) / (t2 - t1));
    Vec3 a3 = p2 * ((t3 - tt) / (t3 - t2)) + p3 * ((tt - t2) / (t3 - t2));
    Vec3 b1 = a1 * ((t2 - tt) / (t2 - t0)) + a2 * ((tt - t0) / (t2 - t0));
    Vec3 b2 = a2 * ((t3 - tt) / (t3 - t1)) + a3 * ((tt - t1) / (t3 - t1));
    return b1 * ((t2 - tt) / (t2 - t1)) + b2 * ((tt - t1) / (t2 - t1));
}

} // namespace bromath
