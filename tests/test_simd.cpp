#include "test_framework.h"

#include "bromath/simd.h"
#include "bromath/mat.h"
#include "bromath/vec.h"

#include <vector>

using namespace bromath;

TEST(simd_float4_basic_ops) {
    Simd4f a(1.0f, 2.0f, 3.0f, 4.0f);
    Simd4f b(5.0f, 6.0f, 7.0f, 8.0f);

    Simd4f sum = a + b;
    ASSERT(nearly(sum.get(0), 6.0f), "simd add 0");
    ASSERT(nearly(sum.get(1), 8.0f), "simd add 1");
    ASSERT(nearly(sum.get(2), 10.0f), "simd add 2");
    ASSERT(nearly(sum.get(3), 12.0f), "simd add 3");

    Simd4f diff = b - a;
    ASSERT(nearly(diff.get(0), 4.0f), "simd sub 0");
    ASSERT(nearly(diff.get(1), 4.0f), "simd sub 1");
    ASSERT(nearly(diff.get(2), 4.0f), "simd sub 2");
    ASSERT(nearly(diff.get(3), 4.0f), "simd sub 3");

    Simd4f prod = a * b;
    ASSERT(nearly(prod.get(0), 5.0f), "simd mul 0");
    ASSERT(nearly(prod.get(1), 12.0f), "simd mul 1");
    ASSERT(nearly(prod.get(2), 21.0f), "simd mul 2");
    ASSERT(nearly(prod.get(3), 32.0f), "simd mul 3");

    Simd4f quot = b / a;
    ASSERT(nearly(quot.get(0), 5.0f), "simd div 0");
    ASSERT(nearly(quot.get(1), 3.0f), "simd div 1");
    ASSERT(nearly(quot.get(2), 7.0f / 3.0f), "simd div 2");
    ASSERT(nearly(quot.get(3), 2.0f), "simd div 3");

    Simd4f fma = simdFmadd(a, b, Simd4f(1.0f));
    ASSERT(nearly(fma.get(0), 6.0f), "simd fma 0");
    ASSERT(nearly(fma.get(1), 13.0f), "simd fma 1");
    ASSERT(nearly(fma.get(2), 22.0f), "simd fma 2");
    ASSERT(nearly(fma.get(3), 33.0f), "simd fma 3");

    Simd4f s = simdSqrt(Simd4f(4.0f, 9.0f, 16.0f, 25.0f));
    ASSERT(nearly(s.get(0), 2.0f), "simd sqrt 0");
    ASSERT(nearly(s.get(1), 3.0f), "simd sqrt 1");
    ASSERT(nearly(s.get(2), 4.0f), "simd sqrt 2");
    ASSERT(nearly(s.get(3), 5.0f), "simd sqrt 3");

    float hsum = simdHadd(a);
    ASSERT(nearly(hsum, 10.0f), "simd hadd");
}

TEST(simd_batch_vector_math) {
    const size_t N = 100;
    std::vector<Vec3> a(N), b(N), out(N);
    for (size_t i = 0; i < N; ++i) {
        a[i] = Vec3{static_cast<float>(i), static_cast<float>(i * 2), static_cast<float>(i * 3)};
        b[i] = Vec3{1.0f, 2.0f, 3.0f};
    }

    batchAdd(a.data(), b.data(), out.data(), N);
    for (size_t i = 0; i < N; ++i) {
        ASSERT(nearly(out[i].x, a[i].x + 1.0f), "batch add x");
        ASSERT(nearly(out[i].y, a[i].y + 2.0f), "batch add y");
        ASSERT(nearly(out[i].z, a[i].z + 3.0f), "batch add z");
    }

    batchScale(a.data(), 0.5f, out.data(), N);
    for (size_t i = 0; i < N; ++i) {
        ASSERT(nearly(out[i].x, a[i].x * 0.5f), "batch scale x");
        ASSERT(nearly(out[i].y, a[i].y * 0.5f), "batch scale y");
        ASSERT(nearly(out[i].z, a[i].z * 0.5f), "batch scale z");
    }

    std::vector<float> dists(N);
    batchDist2ToPoint(a.data(), Vec3{0, 0, 0}, dists.data(), N);
    for (size_t i = 0; i < N; ++i) {
        float expected = vlen2(a[i]);
        ASSERT(nearly(dists[i], expected), "batch dist2 to point");
    }
}

TEST(simd_batch_transform_points_and_vectors) {
    Mat4 m = mtranslate(Vec3{10, 20, 30});
    std::vector<Vec3> pts = {
        {0, 0, 0},
        {1, 2, 3},
        {-5, 10, -15}
    };
    std::vector<Vec3> transformedPts(pts.size());
    std::vector<Vec3> transformedDirs(pts.size());

    batchTransformPoints(m, pts.data(), transformedPts.data(), pts.size());
    for (size_t i = 0; i < pts.size(); ++i) {
        Vec3 expected = pts[i] + Vec3{10, 20, 30};
        ASSERT(nearly(transformedPts[i].x, expected.x), "transformed point x");
        ASSERT(nearly(transformedPts[i].y, expected.y), "transformed point y");
        ASSERT(nearly(transformedPts[i].z, expected.z), "transformed point z");
    }

    batchTransformVectors(m, pts.data(), transformedDirs.data(), pts.size());
    for (size_t i = 0; i < pts.size(); ++i) {
        ASSERT(nearly(transformedDirs[i].x, pts[i].x), "transformed dir x (no translation)");
        ASSERT(nearly(transformedDirs[i].y, pts[i].y), "transformed dir y (no translation)");
        ASSERT(nearly(transformedDirs[i].z, pts[i].z), "transformed dir z (no translation)");
    }
}
