#pragma once

// SIMD and batch vector/matrix math acceleration routines for bromath.
// Automatically leverages SSE2 / AVX2 on x86/x64 and ARM NEON on ARM platforms,
// with clean, robust scalar fallbacks when hardware intrinsics are unavailable.

#include "bromath/aabb.h"
#include "bromath/mat.h"
#include "bromath/scalar.h"
#include "bromath/sphere.h"
#include "bromath/vec.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

#if defined(__AVX2__)
    #include <immintrin.h>
    #define BROMATH_SIMD_AVX2 1
    #define BROMATH_SIMD_SSE2 1
#elif defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
    #include <emmintrin.h>
    #if defined(__SSE4_1__) || defined(_M_X64)
        #include <smmintrin.h>
    #endif
    #define BROMATH_SIMD_SSE2 1
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
    #include <arm_neon.h>
    #define BROMATH_SIMD_NEON 1
#else
    #define BROMATH_SIMD_SCALAR 1
#endif

namespace bromath {

// 4-wide float SIMD vector wrapper.
struct Simd4f {
#if defined(BROMATH_SIMD_SSE2)
    __m128 v;

    Simd4f() = default;
    explicit Simd4f(__m128 val) : v(val) {}
    explicit Simd4f(float s) : v(_mm_set1_ps(s)) {}
    Simd4f(float a, float b, float c, float d) : v(_mm_setr_ps(a, b, c, d)) {}

    static Simd4f load(const float* ptr) {
        return Simd4f(_mm_loadu_ps(ptr));
    }

    void store(float* ptr) const {
        _mm_storeu_ps(ptr, v);
    }

    float get(int idx) const {
        alignas(16) float tmp[4];
        _mm_storeu_ps(tmp, v);
        return tmp[idx & 3];
    }
#elif defined(BROMATH_SIMD_NEON)
    float32x4_t v;

    Simd4f() = default;
    explicit Simd4f(float32x4_t val) : v(val) {}
    explicit Simd4f(float s) : v(vdupq_n_f32(s)) {}
    Simd4f(float a, float b, float c, float d) {
        float tmp[4] = {a, b, c, d};
        v = vld1q_f32(tmp);
    }

    static Simd4f load(const float* ptr) {
        return Simd4f(vld1q_f32(ptr));
    }

    void store(float* ptr) const {
        vst1q_f32(ptr, v);
    }

    float get(int idx) const {
        float tmp[4];
        vst1q_f32(tmp, v);
        return tmp[idx & 3];
    }
#else
    float data[4];

    Simd4f() = default;
    explicit Simd4f(float s) {
        data[0] = s; data[1] = s; data[2] = s; data[3] = s;
    }
    Simd4f(float a, float b, float c, float d) {
        data[0] = a; data[1] = b; data[2] = c; data[3] = d;
    }

    static Simd4f load(const float* ptr) {
        return Simd4f(ptr[0], ptr[1], ptr[2], ptr[3]);
    }

    void store(float* ptr) const {
        ptr[0] = data[0]; ptr[1] = data[1]; ptr[2] = data[2]; ptr[3] = data[3];
    }

    float get(int idx) const {
        return data[idx & 3];
    }
#endif
};

inline Simd4f operator+(Simd4f a, Simd4f b) {
#if defined(BROMATH_SIMD_SSE2)
    return Simd4f(_mm_add_ps(a.v, b.v));
#elif defined(BROMATH_SIMD_NEON)
    return Simd4f(vaddq_f32(a.v, b.v));
#else
    return Simd4f(a.data[0] + b.data[0], a.data[1] + b.data[1],
                  a.data[2] + b.data[2], a.data[3] + b.data[3]);
#endif
}

inline Simd4f operator-(Simd4f a, Simd4f b) {
#if defined(BROMATH_SIMD_SSE2)
    return Simd4f(_mm_sub_ps(a.v, b.v));
#elif defined(BROMATH_SIMD_NEON)
    return Simd4f(vsubq_f32(a.v, b.v));
#else
    return Simd4f(a.data[0] - b.data[0], a.data[1] - b.data[1],
                  a.data[2] - b.data[2], a.data[3] - b.data[3]);
#endif
}

inline Simd4f operator*(Simd4f a, Simd4f b) {
#if defined(BROMATH_SIMD_SSE2)
    return Simd4f(_mm_mul_ps(a.v, b.v));
#elif defined(BROMATH_SIMD_NEON)
    return Simd4f(vmulq_f32(a.v, b.v));
#else
    return Simd4f(a.data[0] * b.data[0], a.data[1] * b.data[1],
                  a.data[2] * b.data[2], a.data[3] * b.data[3]);
#endif
}

inline Simd4f operator/(Simd4f a, Simd4f b) {
#if defined(BROMATH_SIMD_SSE2)
    return Simd4f(_mm_div_ps(a.v, b.v));
#elif defined(BROMATH_SIMD_NEON)
    return Simd4f(vdivq_f32(a.v, b.v));
#else
    return Simd4f(a.data[0] / b.data[0], a.data[1] / b.data[1],
                  a.data[2] / b.data[2], a.data[3] / b.data[3]);
#endif
}

inline Simd4f simdFmadd(Simd4f a, Simd4f b, Simd4f c) {
#if defined(BROMATH_SIMD_AVX2)
    return Simd4f(_mm_fmadd_ps(a.v, b.v, c.v));
#elif defined(BROMATH_SIMD_NEON)
    return Simd4f(vmlaq_f32(c.v, a.v, b.v));
#else
    return (a * b) + c;
#endif
}

inline Simd4f simdMin(Simd4f a, Simd4f b) {
#if defined(BROMATH_SIMD_SSE2)
    return Simd4f(_mm_min_ps(a.v, b.v));
#elif defined(BROMATH_SIMD_NEON)
    return Simd4f(vminq_f32(a.v, b.v));
#else
    return Simd4f(std::min(a.data[0], b.data[0]), std::min(a.data[1], b.data[1]),
                  std::min(a.data[2], b.data[2]), std::min(a.data[3], b.data[3]));
#endif
}

inline Simd4f simdMax(Simd4f a, Simd4f b) {
#if defined(BROMATH_SIMD_SSE2)
    return Simd4f(_mm_max_ps(a.v, b.v));
#elif defined(BROMATH_SIMD_NEON)
    return Simd4f(vmaxq_f32(a.v, b.v));
#else
    return Simd4f(std::max(a.data[0], b.data[0]), std::max(a.data[1], b.data[1]),
                  std::max(a.data[2], b.data[2]), std::max(a.data[3], b.data[3]));
#endif
}

inline Simd4f simdSqrt(Simd4f a) {
#if defined(BROMATH_SIMD_SSE2)
    return Simd4f(_mm_sqrt_ps(a.v));
#elif defined(BROMATH_SIMD_NEON)
    return Simd4f(vsqrtq_f32(a.v));
#else
    return Simd4f(std::sqrt(a.data[0]), std::sqrt(a.data[1]),
                  std::sqrt(a.data[2]), std::sqrt(a.data[3]));
#endif
}

inline float simdHadd(Simd4f a) {
#if defined(BROMATH_SIMD_SSE2)
    __m128 shuf = _mm_shuffle_ps(a.v, a.v, _MM_SHUFFLE(2, 3, 0, 1));
    __m128 sums = _mm_add_ps(a.v, shuf);
    shuf        = _mm_movehl_ps(shuf, sums);
    sums        = _mm_add_ss(sums, shuf);
    return _mm_cvtss_f32(sums);
#elif defined(BROMATH_SIMD_NEON)
    return vaddvq_f32(a.v);
#else
    return a.data[0] + a.data[1] + a.data[2] + a.data[3];
#endif
}

// --- Batch Vector Math Operations ---

inline void batchAdd(const Vec3* a, const Vec3* b, Vec3* out, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        out[i] = a[i] + b[i];
    }
}

inline void batchSub(const Vec3* a, const Vec3* b, Vec3* out, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        out[i] = a[i] - b[i];
    }
}

inline void batchScale(const Vec3* in, float scale, Vec3* out, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        out[i] = in[i] * scale;
    }
}

inline void batchMulAdd(const Vec3* a, float s, const Vec3* b, Vec3* out, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        out[i] = a[i] * s + b[i];
    }
}

inline void batchLerp(const Vec3* a, const Vec3* b, float t, Vec3* out, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        out[i] = vlerp(a[i], b[i], t);
    }
}

inline void batchDot(const Vec3* a, const Vec3* b, float* out, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        out[i] = vdot(a[i], b[i]);
    }
}

inline void batchDist2(const Vec3* a, const Vec3* b, float* out, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        out[i] = vdist2(a[i], b[i]);
    }
}

inline void batchDist2ToPoint(const Vec3* pts, Vec3 target, float* out, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        out[i] = vdist2(pts[i], target);
    }
}

inline void batchDistToPoint(const Vec3* pts, Vec3 target, float* out, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        out[i] = vdist(pts[i], target);
    }
}

inline void batchMin(const Vec3* a, const Vec3* b, Vec3* out, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        out[i] = vmin(a[i], b[i]);
    }
}

inline void batchMax(const Vec3* a, const Vec3* b, Vec3* out, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        out[i] = vmax(a[i], b[i]);
    }
}

inline void batchNormalize(const Vec3* in, Vec3* out, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        out[i] = vnorm(in[i]);
    }
}

// Transform a batch of 3D points (assumes w=1) using a 4x4 matrix.
inline void batchTransformPoints(const Mat4& m, const Vec3* in, Vec3* out, size_t count) {
#if defined(BROMATH_SIMD_SSE2)
    __m128 c0 = _mm_setr_ps(m.at(0, 0), m.at(1, 0), m.at(2, 0), 0.0f);
    __m128 c1 = _mm_setr_ps(m.at(0, 1), m.at(1, 1), m.at(2, 1), 0.0f);
    __m128 c2 = _mm_setr_ps(m.at(0, 2), m.at(1, 2), m.at(2, 2), 0.0f);
    __m128 c3 = _mm_setr_ps(m.at(0, 3), m.at(1, 3), m.at(2, 3), 0.0f);

    for (size_t i = 0; i < count; ++i) {
        __m128 x = _mm_set1_ps(in[i].x);
        __m128 y = _mm_set1_ps(in[i].y);
        __m128 z = _mm_set1_ps(in[i].z);

        __m128 res = _mm_add_ps(_mm_mul_ps(x, c0), _mm_mul_ps(y, c1));
        res = _mm_add_ps(res, _mm_mul_ps(z, c2));
        res = _mm_add_ps(res, c3);

        alignas(16) float tmp[4];
        _mm_storeu_ps(tmp, res);
        out[i] = Vec3{tmp[0], tmp[1], tmp[2]};
    }
#else
    for (size_t i = 0; i < count; ++i) {
        out[i] = mtransformPoint(m, in[i]);
    }
#endif
}

// Transform a batch of 3D direction vectors (assumes w=0) using a 4x4 matrix.
inline void batchTransformVectors(const Mat4& m, const Vec3* in, Vec3* out, size_t count) {
#if defined(BROMATH_SIMD_SSE2)
    __m128 c0 = _mm_setr_ps(m.at(0, 0), m.at(1, 0), m.at(2, 0), 0.0f);
    __m128 c1 = _mm_setr_ps(m.at(0, 1), m.at(1, 1), m.at(2, 1), 0.0f);
    __m128 c2 = _mm_setr_ps(m.at(0, 2), m.at(1, 2), m.at(2, 2), 0.0f);

    for (size_t i = 0; i < count; ++i) {
        __m128 x = _mm_set1_ps(in[i].x);
        __m128 y = _mm_set1_ps(in[i].y);
        __m128 z = _mm_set1_ps(in[i].z);

        __m128 res = _mm_add_ps(_mm_mul_ps(x, c0), _mm_mul_ps(y, c1));
        res = _mm_add_ps(res, _mm_mul_ps(z, c2));

        alignas(16) float tmp[4];
        _mm_storeu_ps(tmp, res);
        out[i] = Vec3{tmp[0], tmp[1], tmp[2]};
    }
#else
    for (size_t i = 0; i < count; ++i) {
        out[i] = mtransformDir(m, in[i]);
    }
#endif
}

// Batch test whether spheres overlap a target sphere (1 if overlap, 0 if not).
inline void batchSphereOverlap(const Sphere* spheres, Sphere target, uint8_t* outMask, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        float rSum = spheres[i].radius + target.radius;
        outMask[i] = (vdist2(spheres[i].center, target.center) <= rSum * rSum) ? 1 : 0;
    }
}

// Batch test whether points lie inside an AABB (1 if inside, 0 if not).
inline void batchAABBContains(const AABB3& box, const Vec3* pts, uint8_t* outMask, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        outMask[i] = acontains(box, pts[i]) ? 1 : 0;
    }
}

// Find all points within radius of center, appending matching indices.
inline void simdFindWithinRadius(const Vec3* pts, size_t count, Vec3 center, float radius, std::vector<uint32_t>& outIndices) {
    float r2 = radius * radius;
    for (size_t i = 0; i < count; ++i) {
        if (vdist2(pts[i], center) <= r2) {
            outIndices.push_back(static_cast<uint32_t>(i));
        }
    }
}

} // namespace bromath
