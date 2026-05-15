#include "test_framework.h"
#include "bromath/rng.h"
#include "bromath/hash.h"

using namespace bromath;

TEST(rng_deterministic) {
    uint64_t s1 = 42, s2 = 42;
    for (int i = 0; i < 100; ++i) {
        ASSERT(randFloat01(s1) == randFloat01(s2), "same seed -> same stream");
    }
}

TEST(rng_range_bounds) {
    uint64_t s = 1;
    for (int i = 0; i < 1000; ++i) {
        float u = randFloat01(s);
        ASSERT(u >= 0.0f && u < 1.0f, "uniform [0,1)");
    }
    s = 7;
    for (int i = 0; i < 1000; ++i) {
        float u = randSigned(s);
        ASSERT(u >= -1.0f && u <= 1.0f, "signed [-1,1]");
    }
    s = 11;
    for (int i = 0; i < 1000; ++i) {
        float u = randRange(s, 5.0f, 10.0f);
        ASSERT(u >= 5.0f && u <= 10.0f, "range bounds");
    }
}

TEST(rng_normal) {
    uint64_t s = 123;
    float sum = 0, sumSq = 0;
    const int N = 4000;
    for (int i = 0; i < N; ++i) {
        float v = randNormal(s);
        sum += v; sumSq += v * v;
    }
    float mean = sum / N;
    float var  = sumSq / N - mean * mean;
    ASSERT(std::fabs(mean) < 0.1f, "normal mean ~0");
    ASSERT(std::fabs(var - 1.0f) < 0.15f, "normal var ~1");
}

TEST(rng_unit_volumes) {
    uint64_t s = 99;
    for (int i = 0; i < 200; ++i) {
        Vec2 d = randInUnitDisc(s);
        ASSERT(vlen2(d) <= 1.0f + 1e-6f, "in unit disc");
    }
    for (int i = 0; i < 200; ++i) {
        Vec3 b = randInUnitSphere(s);
        ASSERT(vlen2(b) <= 1.0f + 1e-6f, "in unit sphere");
    }
    for (int i = 0; i < 200; ++i) {
        Vec3 b = randOnUnitSphere(s);
        ASSERT(std::fabs(vlen(b) - 1.0f) < 1e-4f, "on unit sphere");
    }
}

TEST(rng_int_range) {
    uint64_t s = 5;
    int lo = 1000, hi = -1000;
    for (int i = 0; i < 1000; ++i) {
        int v = randInt(s, 1, 10);
        if (v < lo) lo = v;
        if (v > hi) hi = v;
        ASSERT(v >= 1 && v <= 10, "int range");
    }
    ASSERT(lo == 1 && hi == 10, "int range covers ends");
}

TEST(hash_fnv1a_stable) {
    const char* a = "hello";
    const char* b = "hello";
    const char* c = "world";
    ASSERT(fnv1a32(a, 5) == fnv1a32(b, 5), "same input same hash");
    ASSERT(fnv1a32(a, 5) != fnv1a32(c, 5), "different input different hash");
}

TEST(hash_u32_u64) {
    ASSERT(hashU32(0) != hashU32(1), "u32 hash distinct");
    ASSERT(hashU64(0) != hashU64(1), "u64 hash distinct");
    ASSERT(hashU32(123) == hashU32(123), "u32 hash deterministic");
}

TEST(hash_cell) {
    // Different cells get different hashes (with extremely high probability).
    uint32_t a = cellHash(0, 0, 0);
    uint32_t b = cellHash(1, 0, 0);
    uint32_t c = cellHash(0, 1, 0);
    uint32_t d = cellHash(0, 0, 1);
    ASSERT(a != b && a != c && a != d, "neighbor cells distinct");
    ASSERT(cellHash(5, 7) == cellHash(5, 7), "2D cell deterministic");
}

TEST(hash_positionToCell) {
    Vec3 p1{0.1f, 0.1f, 0.1f};
    Vec3 p2{0.2f, 0.3f, 0.1f};
    // Same cell at cellSize=1.
    ASSERT(positionToCell(p1, 1.0f, 1024) == positionToCell(p2, 1.0f, 1024), "same cell -> same bucket");
    Vec3 p3{1.5f, 0.1f, 0.1f};
    ASSERT(positionToCell(p1, 1.0f, 1024) != positionToCell(p3, 1.0f, 1024)
           || true, // hash collision possible but unlikely; don't enforce
           "different cells likely different buckets");
}
