#pragma once

// View frustum: six bounding planes (left, right, bottom, top, near, far)
// with inward-facing normals. A point is inside the frustum iff every plane
// returns a non-negative signed distance.
//
// Extracted from a view-projection matrix via Gribb-Hartmann: each plane is
// a linear combination of matrix rows. Works for both GL and DX style VP
// matrices since we normalize the planes.

#include "bromath/aabb.h"
#include "bromath/mat.h"
#include "bromath/plane.h"
#include "bromath/sphere.h"
#include "bromath/vec.h"

#include <cmath>

namespace bromath {

struct Frustum {
    Plane planes[6]; // left, right, bottom, top, near, far
};

inline Frustum ffromViewProj(const Mat4& vp) {
    auto row = [&](int r) {
        return Vec3{ vp.at(r, 0), vp.at(r, 1), vp.at(r, 2) };
    };
    auto row_w = [&](int r) { return vp.at(r, 3); };

    auto makePlane = [&](Vec3 n, float d) {
        float L = vlen(n);
        if (L < 1e-20f) return Plane{ {0, 1, 0}, 0 };
        float inv = 1.0f / L;
        return Plane{ n * inv, d * inv };
    };

    Vec3 r0 = row(0), r1 = row(1), r2 = row(2), r3 = row(3);
    float w0 = row_w(0), w1 = row_w(1), w2 = row_w(2), w3 = row_w(3);

    Frustum f;
    f.planes[0] = makePlane(r3 + r0, w3 + w0); // left
    f.planes[1] = makePlane(r3 - r0, w3 - w0); // right
    f.planes[2] = makePlane(r3 + r1, w3 + w1); // bottom
    f.planes[3] = makePlane(r3 - r1, w3 - w1); // top
    f.planes[4] = makePlane(r3 + r2, w3 + w2); // near
    f.planes[5] = makePlane(r3 - r2, w3 - w2); // far
    return f;
}

inline bool fcontains(const Frustum& f, Vec3 p) {
    for (int i = 0; i < 6; ++i) {
        if (psignedDistance(f.planes[i], p) < 0.0f) return false;
    }
    return true;
}

// Conservative AABB cull: returns false only if the box is fully outside
// at least one plane. May return true for boxes that straddle the frustum
// corners (standard "false positives" for plane-by-plane culling).
inline bool fintersects(const Frustum& f, const AABB3& a) {
    for (int i = 0; i < 6; ++i) {
        const Plane& p = f.planes[i];
        // The "positive vertex" of the AABB relative to the plane normal —
        // pick the corner farthest along the normal. If it's behind the
        // plane, the entire box is behind.
        Vec3 pv{
            p.normal.x >= 0.0f ? a.max.x : a.min.x,
            p.normal.y >= 0.0f ? a.max.y : a.min.y,
            p.normal.z >= 0.0f ? a.max.z : a.min.z
        };
        if (psignedDistance(p, pv) < 0.0f) return false;
    }
    return true;
}

inline bool fintersects(const Frustum& f, const Sphere& s) {
    for (int i = 0; i < 6; ++i) {
        if (psignedDistance(f.planes[i], s.center) < -s.radius) return false;
    }
    return true;
}

} // namespace bromath
