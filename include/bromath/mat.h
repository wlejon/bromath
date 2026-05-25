#pragma once

// Mat4: 4x4 column-major (OpenGL / glTF / Jolt convention). Stored as
// 16 floats with the column index varying slowest:
//
//     m.data[col * 4 + row]
//
// equivalent to:
//
//     | m00 m01 m02 m03 |    where m(row, col) = data[col * 4 + row]
//     | m10 m11 m12 m13 |
//     | m20 m21 m22 m23 |
//     | m30 m31 m32 m33 |
//
// data() returns a pointer suitable for glUniformMatrix4fv (transpose=GL_FALSE).

#include "bromath/quat.h"
#include "bromath/scalar.h"
#include "bromath/vec.h"

#include <cmath>

namespace bromath {

struct Mat4 {
    float data[16] = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    };

    constexpr float& at(int row, int col)       { return data[col * 4 + row]; }
    constexpr float  at(int row, int col) const { return data[col * 4 + row]; }
};

inline constexpr Mat4 midentity() { return Mat4{}; }

inline constexpr Mat4 mmul(const Mat4& a, const Mat4& b) {
    Mat4 r;
    for (int c = 0; c < 4; ++c) {
        for (int rr = 0; rr < 4; ++rr) {
            float s = 0.0f;
            for (int k = 0; k < 4; ++k) {
                s += a.at(rr, k) * b.at(k, c);
            }
            r.at(rr, c) = s;
        }
    }
    return r;
}

inline constexpr Mat4 mtranspose(const Mat4& m) {
    Mat4 r;
    for (int c = 0; c < 4; ++c)
        for (int rr = 0; rr < 4; ++rr)
            r.at(rr, c) = m.at(c, rr);
    return r;
}

// General 4x4 inverse via cofactor expansion. Returns identity if singular.
inline Mat4 minverse(const Mat4& m) {
    const float* a = m.data;
    float inv[16];
    inv[0]  =  a[5]*a[10]*a[15] - a[5]*a[11]*a[14] - a[9]*a[6]*a[15]
             + a[9]*a[7]*a[14]  + a[13]*a[6]*a[11] - a[13]*a[7]*a[10];
    inv[4]  = -a[4]*a[10]*a[15] + a[4]*a[11]*a[14] + a[8]*a[6]*a[15]
             - a[8]*a[7]*a[14]  - a[12]*a[6]*a[11] + a[12]*a[7]*a[10];
    inv[8]  =  a[4]*a[9]*a[15]  - a[4]*a[11]*a[13] - a[8]*a[5]*a[15]
             + a[8]*a[7]*a[13]  + a[12]*a[5]*a[11] - a[12]*a[7]*a[9];
    inv[12] = -a[4]*a[9]*a[14]  + a[4]*a[10]*a[13] + a[8]*a[5]*a[14]
             - a[8]*a[6]*a[13]  - a[12]*a[5]*a[10] + a[12]*a[6]*a[9];
    inv[1]  = -a[1]*a[10]*a[15] + a[1]*a[11]*a[14] + a[9]*a[2]*a[15]
             - a[9]*a[3]*a[14]  - a[13]*a[2]*a[11] + a[13]*a[3]*a[10];
    inv[5]  =  a[0]*a[10]*a[15] - a[0]*a[11]*a[14] - a[8]*a[2]*a[15]
             + a[8]*a[3]*a[14]  + a[12]*a[2]*a[11] - a[12]*a[3]*a[10];
    inv[9]  = -a[0]*a[9]*a[15]  + a[0]*a[11]*a[13] + a[8]*a[1]*a[15]
             - a[8]*a[3]*a[13]  - a[12]*a[1]*a[11] + a[12]*a[3]*a[9];
    inv[13] =  a[0]*a[9]*a[14]  - a[0]*a[10]*a[13] - a[8]*a[1]*a[14]
             + a[8]*a[2]*a[13]  + a[12]*a[1]*a[10] - a[12]*a[2]*a[9];
    inv[2]  =  a[1]*a[6]*a[15]  - a[1]*a[7]*a[14]  - a[5]*a[2]*a[15]
             + a[5]*a[3]*a[14]  + a[13]*a[2]*a[7]  - a[13]*a[3]*a[6];
    inv[6]  = -a[0]*a[6]*a[15]  + a[0]*a[7]*a[14]  + a[4]*a[2]*a[15]
             - a[4]*a[3]*a[14]  - a[12]*a[2]*a[7]  + a[12]*a[3]*a[6];
    inv[10] =  a[0]*a[5]*a[15]  - a[0]*a[7]*a[13]  - a[4]*a[1]*a[15]
             + a[4]*a[3]*a[13]  + a[12]*a[1]*a[7]  - a[12]*a[3]*a[5];
    inv[14] = -a[0]*a[5]*a[14]  + a[0]*a[6]*a[13]  + a[4]*a[1]*a[14]
             - a[4]*a[2]*a[13]  - a[12]*a[1]*a[6]  + a[12]*a[2]*a[5];
    inv[3]  = -a[1]*a[6]*a[11]  + a[1]*a[7]*a[10]  + a[5]*a[2]*a[11]
             - a[5]*a[3]*a[10]  - a[9]*a[2]*a[7]   + a[9]*a[3]*a[6];
    inv[7]  =  a[0]*a[6]*a[11]  - a[0]*a[7]*a[10]  - a[4]*a[2]*a[11]
             + a[4]*a[3]*a[10]  + a[8]*a[2]*a[7]   - a[8]*a[3]*a[6];
    inv[11] = -a[0]*a[5]*a[11]  + a[0]*a[7]*a[9]   + a[4]*a[1]*a[11]
             - a[4]*a[3]*a[9]   - a[8]*a[1]*a[7]   + a[8]*a[3]*a[5];
    inv[15] =  a[0]*a[5]*a[10]  - a[0]*a[6]*a[9]   - a[4]*a[1]*a[10]
             + a[4]*a[2]*a[9]   + a[8]*a[1]*a[6]   - a[8]*a[2]*a[5];

    float det = a[0]*inv[0] + a[1]*inv[4] + a[2]*inv[8] + a[3]*inv[12];
    if (std::fabs(det) < 1e-20f) return midentity();
    float invDet = 1.0f / det;
    Mat4 r;
    for (int i = 0; i < 16; ++i) r.data[i] = inv[i] * invDet;
    return r;
}

inline constexpr Mat4 mtranslate(Vec3 t) {
    Mat4 m;
    m.at(0, 3) = t.x;
    m.at(1, 3) = t.y;
    m.at(2, 3) = t.z;
    return m;
}

inline constexpr Mat4 mscale(Vec3 s) {
    Mat4 m;
    m.at(0, 0) = s.x;
    m.at(1, 1) = s.y;
    m.at(2, 2) = s.z;
    return m;
}

// Build a rotation Mat4 from a unit quaternion.
inline constexpr Mat4 mfromQuat(Quat q) {
    float xx = q.x*q.x, yy = q.y*q.y, zz = q.z*q.z;
    float xy = q.x*q.y, xz = q.x*q.z, yz = q.y*q.z;
    float wx = q.w*q.x, wy = q.w*q.y, wz = q.w*q.z;
    Mat4 m;
    m.at(0,0) = 1 - 2*(yy + zz); m.at(0,1) = 2*(xy - wz);     m.at(0,2) = 2*(xz + wy);
    m.at(1,0) = 2*(xy + wz);     m.at(1,1) = 1 - 2*(xx + zz); m.at(1,2) = 2*(yz - wx);
    m.at(2,0) = 2*(xz - wy);     m.at(2,1) = 2*(yz + wx);     m.at(2,2) = 1 - 2*(xx + yy);
    m.at(3,3) = 1;
    return m;
}

// Compose translation * rotation * scale.
inline constexpr Mat4 mfromTRS(Vec3 t, Quat r, Vec3 s) {
    Mat4 m = mfromQuat(r);
    m.at(0,0) *= s.x; m.at(1,0) *= s.x; m.at(2,0) *= s.x;
    m.at(0,1) *= s.y; m.at(1,1) *= s.y; m.at(2,1) *= s.y;
    m.at(0,2) *= s.z; m.at(1,2) *= s.z; m.at(2,2) *= s.z;
    m.at(0,3) = t.x;
    m.at(1,3) = t.y;
    m.at(2,3) = t.z;
    return m;
}

// Decompose an affine TRS matrix back into components. Skew is discarded.
inline void mdecompose(const Mat4& m, Vec3& t, Quat& r, Vec3& s) {
    t = { m.at(0,3), m.at(1,3), m.at(2,3) };
    Vec3 c0{ m.at(0,0), m.at(1,0), m.at(2,0) };
    Vec3 c1{ m.at(0,1), m.at(1,1), m.at(2,1) };
    Vec3 c2{ m.at(0,2), m.at(1,2), m.at(2,2) };
    s = { vlen(c0), vlen(c1), vlen(c2) };
    if (s.x > 1e-20f) c0 = c0 / s.x;
    if (s.y > 1e-20f) c1 = c1 / s.y;
    if (s.z > 1e-20f) c2 = c2 / s.z;
    // Build a quaternion from the orthonormal basis (Shepperd's method).
    float tr = c0.x + c1.y + c2.z;
    if (tr > 0.0f) {
        float ss = std::sqrt(tr + 1.0f) * 2.0f;
        r.w = 0.25f * ss;
        r.x = (c1.z - c2.y) / ss;
        r.y = (c2.x - c0.z) / ss;
        r.z = (c0.y - c1.x) / ss;
    } else if (c0.x > c1.y && c0.x > c2.z) {
        float ss = std::sqrt(1.0f + c0.x - c1.y - c2.z) * 2.0f;
        r.w = (c1.z - c2.y) / ss;
        r.x = 0.25f * ss;
        r.y = (c1.x + c0.y) / ss;
        r.z = (c2.x + c0.z) / ss;
    } else if (c1.y > c2.z) {
        float ss = std::sqrt(1.0f + c1.y - c0.x - c2.z) * 2.0f;
        r.w = (c2.x - c0.z) / ss;
        r.x = (c1.x + c0.y) / ss;
        r.y = 0.25f * ss;
        r.z = (c2.y + c1.z) / ss;
    } else {
        float ss = std::sqrt(1.0f + c2.z - c0.x - c1.y) * 2.0f;
        r.w = (c0.y - c1.x) / ss;
        r.x = (c2.x + c0.z) / ss;
        r.y = (c2.y + c1.z) / ss;
        r.z = 0.25f * ss;
    }
}

// Right-handed view matrix looking from `eye` toward `center`, with `up`
// approximating world up. Matches gluLookAt.
inline Mat4 mlookAt(Vec3 eye, Vec3 center, Vec3 up) {
    Vec3 f = vnorm(center - eye);          // forward
    Vec3 s = vnorm(vcross(f, up));         // right
    Vec3 u = vcross(s, f);                 // recomputed up
    Mat4 m;
    m.at(0,0) = s.x;  m.at(0,1) = s.y;  m.at(0,2) = s.z;  m.at(0,3) = -vdot(s, eye);
    m.at(1,0) = u.x;  m.at(1,1) = u.y;  m.at(1,2) = u.z;  m.at(1,3) = -vdot(u, eye);
    m.at(2,0) = -f.x; m.at(2,1) = -f.y; m.at(2,2) = -f.z; m.at(2,3) =  vdot(f, eye);
    m.at(3,3) = 1.0f;
    return m;
}

// Right-handed perspective. fovY in radians. Maps NDC z to [-1, 1] (GL).
inline Mat4 mperspective(float fovY, float aspect, float znear, float zfar) {
    float f = 1.0f / std::tan(fovY * 0.5f);
    Mat4 m;
    for (int i = 0; i < 16; ++i) m.data[i] = 0.0f;
    m.at(0,0) = f / aspect;
    m.at(1,1) = f;
    m.at(2,2) = (zfar + znear) / (znear - zfar);
    m.at(2,3) = (2.0f * zfar * znear) / (znear - zfar);
    m.at(3,2) = -1.0f;
    return m;
}

// Right-handed orthographic. NDC z in [-1, 1] (GL).
inline Mat4 mortho(float l, float r, float b, float t, float znear, float zfar) {
    Mat4 m;
    for (int i = 0; i < 16; ++i) m.data[i] = 0.0f;
    m.at(0,0) =  2.0f / (r - l);
    m.at(1,1) =  2.0f / (t - b);
    m.at(2,2) = -2.0f / (zfar - znear);
    m.at(0,3) = -(r + l) / (r - l);
    m.at(1,3) = -(t + b) / (t - b);
    m.at(2,3) = -(zfar + znear) / (zfar - znear);
    m.at(3,3) = 1.0f;
    return m;
}

// Transform a point (assumes w=1). Includes translation.
inline constexpr Vec3 mtransformPoint(const Mat4& m, Vec3 p) {
    return {
        m.at(0,0)*p.x + m.at(0,1)*p.y + m.at(0,2)*p.z + m.at(0,3),
        m.at(1,0)*p.x + m.at(1,1)*p.y + m.at(1,2)*p.z + m.at(1,3),
        m.at(2,0)*p.x + m.at(2,1)*p.y + m.at(2,2)*p.z + m.at(2,3)
    };
}

// Transform a direction (assumes w=0). Skips translation.
inline constexpr Vec3 mtransformDir(const Mat4& m, Vec3 d) {
    return {
        m.at(0,0)*d.x + m.at(0,1)*d.y + m.at(0,2)*d.z,
        m.at(1,0)*d.x + m.at(1,1)*d.y + m.at(1,2)*d.z,
        m.at(2,0)*d.x + m.at(2,1)*d.y + m.at(2,2)*d.z
    };
}

} // namespace bromath
