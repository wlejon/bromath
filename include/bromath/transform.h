#pragma once

// Transform: position + rotation + scale. The natural representation for
// scene nodes; converts to/from Mat4 on demand.

#include "bromath/mat.h"
#include "bromath/quat.h"
#include "bromath/vec.h"

namespace bromath {

struct Transform {
    Vec3 position = {0, 0, 0};
    Quat rotation = {0, 0, 0, 1};
    Vec3 scale    = {1, 1, 1};
};

inline constexpr Transform tidentity() { return {}; }

inline Mat4 ttoMat4(const Transform& t) {
    return mfromTRS(t.position, t.rotation, t.scale);
}

inline Transform tfromMat4(const Mat4& m) {
    Transform t;
    mdecompose(m, t.position, t.rotation, t.scale);
    return t;
}

// Compose two transforms: result = parent then child (parent applied first).
inline Transform tmul(const Transform& parent, const Transform& child) {
    Transform r;
    r.scale    = { parent.scale.x * child.scale.x,
                   parent.scale.y * child.scale.y,
                   parent.scale.z * child.scale.z };
    r.rotation = qmul(parent.rotation, child.rotation);
    Vec3 scaledChildPos = {
        child.position.x * parent.scale.x,
        child.position.y * parent.scale.y,
        child.position.z * parent.scale.z
    };
    r.position = parent.position + qrotate(parent.rotation, scaledChildPos);
    return r;
}

inline Vec3 ttransformPoint(const Transform& t, Vec3 p) {
    Vec3 scaled = { p.x * t.scale.x, p.y * t.scale.y, p.z * t.scale.z };
    return t.position + qrotate(t.rotation, scaled);
}

inline Vec3 ttransformDir(const Transform& t, Vec3 d) {
    return qrotate(t.rotation, d);
}

} // namespace bromath
