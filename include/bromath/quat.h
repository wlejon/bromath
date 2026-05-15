#pragma once

// Quaternion (xyzw order, matching glTF / Jolt / glm). Identity = (0,0,0,1).
// Free-function ops following the same style as vec.h.

#include "bromath/scalar.h"
#include "bromath/vec.h"

#include <cmath>

namespace bromath {

struct Quat {
    float x = 0, y = 0, z = 0, w = 1;
    constexpr Quat() = default;
    constexpr Quat(float xx, float yy, float zz, float ww) : x(xx), y(yy), z(zz), w(ww) {}
};

inline constexpr Quat qidentity() { return {0, 0, 0, 1}; }

inline constexpr float qdot(Quat a, Quat b) {
    return a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w;
}

inline constexpr float qlen2(Quat q) { return qdot(q, q); }
inline float           qlen(Quat q)  { return std::sqrt(qlen2(q)); }

inline Quat qnorm(Quat q) {
    float L = qlen(q);
    if (L < 1e-20f) return qidentity();
    float inv = 1.0f / L;
    return {q.x*inv, q.y*inv, q.z*inv, q.w*inv};
}

inline constexpr Quat qconjugate(Quat q) { return {-q.x, -q.y, -q.z, q.w}; }
// Inverse of a unit quaternion is its conjugate.
inline constexpr Quat qinverse(Quat q)   { return qconjugate(q); }

inline constexpr Quat qmul(Quat a, Quat b) {
    return {
        a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
        a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
        a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w,
        a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z
    };
}

// Rotate vector `v` by unit quaternion `q`.
inline constexpr Vec3 qrotate(Quat q, Vec3 v) {
    Vec3 u{q.x, q.y, q.z};
    float s = q.w;
    return u * (2.0f * vdot(u, v))
         + v * (s*s - vdot(u, u))
         + vcross(u, v) * (2.0f * s);
}

inline Quat qaxisAngle(Vec3 axis, float angleRad) {
    Vec3 a = vnorm(axis);
    float h = angleRad * 0.5f;
    float s = std::sin(h);
    return {a.x * s, a.y * s, a.z * s, std::cos(h)};
}

// Shortest-arc quaternion that rotates unit vector `from` to unit vector `to`.
inline Quat qfromTo(Vec3 from, Vec3 to) {
    float d = vdot(from, to);
    if (d > 0.999999f) return qidentity();
    if (d < -0.999999f) {
        // 180° flip: pick any axis perpendicular to `from`.
        Vec3 axis = vcross({1, 0, 0}, from);
        if (vlen2(axis) < 1e-12f) axis = vcross({0, 1, 0}, from);
        axis = vnorm(axis);
        return {axis.x, axis.y, axis.z, 0};
    }
    Vec3 c = vcross(from, to);
    float s = std::sqrt((1.0f + d) * 2.0f);
    float invs = 1.0f / s;
    return {c.x * invs, c.y * invs, c.z * invs, s * 0.5f};
}

// Euler-XYZ: roll (x) -> pitch (y) -> yaw (z), intrinsic. Radians.
inline Quat qfromEuler(float rx, float ry, float rz) {
    float cx = std::cos(rx * 0.5f), sx = std::sin(rx * 0.5f);
    float cy = std::cos(ry * 0.5f), sy = std::sin(ry * 0.5f);
    float cz = std::cos(rz * 0.5f), sz = std::sin(rz * 0.5f);
    return {
        sx*cy*cz - cx*sy*sz,
        cx*sy*cz + sx*cy*sz,
        cx*cy*sz - sx*sy*cz,
        cx*cy*cz + sx*sy*sz
    };
}

inline Quat qfromEuler(Vec3 r) { return qfromEuler(r.x, r.y, r.z); }

// Decompose unit quaternion to Euler XYZ (rx, ry, rz). Handles gimbal pole.
inline Vec3 qtoEuler(Quat q) {
    Vec3 e;
    float sinp = 2.0f * (q.w * q.y - q.z * q.x);
    if (std::fabs(sinp) >= 1.0f) {
        // Gimbal lock — clamp pitch, fold roll into yaw.
        e.y = std::copysign(HALF_PI, sinp);
        e.x = 0.0f;
        e.z = std::atan2(-2.0f * (q.x * q.y - q.w * q.z),
                          1.0f - 2.0f * (q.y * q.y + q.z * q.z));
    } else {
        e.y = std::asin(sinp);
        e.x = std::atan2(2.0f * (q.w * q.x + q.y * q.z),
                          1.0f - 2.0f * (q.x * q.x + q.y * q.y));
        e.z = std::atan2(2.0f * (q.w * q.z + q.x * q.y),
                          1.0f - 2.0f * (q.y * q.y + q.z * q.z));
    }
    return e;
}

// Normalized lerp — cheap, not constant-velocity. Fine for short steps.
inline Quat qnlerp(Quat a, Quat b, float t) {
    if (qdot(a, b) < 0.0f) b = {-b.x, -b.y, -b.z, -b.w};
    Quat r{
        lerp(a.x, b.x, t),
        lerp(a.y, b.y, t),
        lerp(a.z, b.z, t),
        lerp(a.w, b.w, t)
    };
    return qnorm(r);
}

// Spherical linear interpolation — constant angular velocity. Falls back to
// nlerp when the quaternions are nearly parallel.
inline Quat qslerp(Quat a, Quat b, float t) {
    float d = qdot(a, b);
    if (d < 0.0f) { b = {-b.x, -b.y, -b.z, -b.w}; d = -d; }
    if (d > 0.9995f) return qnlerp(a, b, t);
    float theta0 = std::acos(d);
    float theta  = theta0 * t;
    float sinT0  = std::sin(theta0);
    float s0 = std::sin(theta0 - theta) / sinT0;
    float s1 = std::sin(theta) / sinT0;
    return {
        a.x*s0 + b.x*s1,
        a.y*s0 + b.y*s1,
        a.z*s0 + b.z*s1,
        a.w*s0 + b.w*s1
    };
}

} // namespace bromath
