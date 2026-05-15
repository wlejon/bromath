#include "test_framework.h"

#include "bromath/aabb.h"
#include "bromath/spatial_hash.h"
#include "bromath/sphere.h"
#include "bromath/vec.h"

#include <algorithm>
#include <random>
#include <set>
#include <vector>

using namespace bromath;

TEST(spatial_hash_radius_query_points) {
    SpatialHash3D hash(1.0f);
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> uni(-10.0f, 10.0f);
    std::vector<Vec3> pts;
    pts.reserve(1000);
    for (int i = 0; i < 1000; ++i) {
        Vec3 p{uni(rng), uni(rng), uni(rng)};
        pts.push_back(p);
        hash.insert(p, i);
    }
    ASSERT(hash.size() == 1000, "size after insert");

    Vec3 center{0, 0, 0};
    float radius = 3.0f;
    std::set<int32_t> brute;
    for (size_t i = 0; i < pts.size(); ++i) {
        if (vdist(pts[i], center) <= radius) brute.insert(static_cast<int32_t>(i));
    }
    std::vector<int32_t> got;
    hash.radiusQuery(center, radius, got);
    std::set<int32_t> gotSet(got.begin(), got.end());
    ASSERT(brute == gotSet, "radius query matches brute force");

    int32_t bestId = -1;
    float bestD2 = 1e30f;
    for (size_t i = 0; i < pts.size(); ++i) {
        float d2 = vdist2(pts[i], center);
        if (d2 < bestD2) { bestD2 = d2; bestId = static_cast<int32_t>(i); }
    }
    int32_t got2 = hash.nearest(center, 100.0f);
    ASSERT(got2 == bestId, "nearest matches brute force");
}

TEST(spatial_hash_remove_clear) {
    SpatialHash3D h(0.5f);
    h.insert(Vec3{0, 0, 0}, 1);
    h.insert(Vec3{0.1f, 0, 0}, 2);
    h.insert(Vec3{5, 5, 5}, 3);
    ASSERT(h.size() == 3, "three inserted");
    h.remove(2);
    ASSERT(h.size() == 2, "size after remove");
    std::vector<int32_t> ids;
    h.radiusQuery(Vec3{0, 0, 0}, 1.0f, ids);
    ASSERT(ids.size() == 1 && ids[0] == 1, "removed id is gone");
    h.clear();
    ASSERT(h.size() == 0, "size after clear");
}

TEST(spatial_hash_sphere_insert_finds_distant_center) {
    // Sphere whose center sits several cells away but whose surface reaches
    // the query origin must be found.
    SpatialHash3D h(0.5f);
    h.insert(Sphere{{5.0f, 0, 0}, 2.0f}, 100); // surface at x=3
    h.insert(Sphere{{5.0f, 0, 0}, 10.0f}, 101); // surface at x=-5
    h.insert(Sphere{{50.0f, 0, 0}, 1.0f}, 102); // far away, should miss

    std::vector<int32_t> hits;
    h.radiusQuery(Vec3{0, 0, 0}, 0.5f, hits);
    std::set<int32_t> set(hits.begin(), hits.end());
    ASSERT(set.count(101) == 1, "huge sphere matches query at origin");
    ASSERT(set.count(100) == 0, "small sphere at x=5 does not reach origin");
    ASSERT(set.count(102) == 0, "far sphere does not reach origin");
    ASSERT(h.maxRadius() == 10.0f, "maxRadius tracked");
}

TEST(spatial_hash_query_aabb_points) {
    SpatialHash3D h(1.0f);
    h.insert(Vec3{0.5f, 0.5f, 0.5f}, 1);
    h.insert(Vec3{1.5f, 0.5f, 0.5f}, 2);
    h.insert(Vec3{10.0f, 10.0f, 10.0f}, 3);

    AABB3 box{{0, 0, 0}, {1, 1, 1}};
    std::vector<int32_t> hits;
    h.queryAABB(box, hits);
    std::set<int32_t> set(hits.begin(), hits.end());
    ASSERT(set.size() == 1 && set.count(1) == 1, "only point inside box is returned");
}

TEST(spatial_hash_query_aabb_spheres) {
    SpatialHash3D h(0.5f);
    // Sphere whose center is outside the box but whose body intrudes.
    h.insert(Sphere{{2.5f, 0.5f, 0.5f}, 2.0f}, 1); // reaches x=0.5, inside
    h.insert(Sphere{{2.5f, 0.5f, 0.5f}, 0.5f}, 2); // reaches x=2.0, outside
    h.insert(Sphere{{0.5f, 0.5f, 0.5f}, 0.1f}, 3); // fully inside

    AABB3 box{{0, 0, 0}, {1, 1, 1}};
    std::vector<int32_t> hits;
    h.queryAABB(box, hits);
    std::set<int32_t> set(hits.begin(), hits.end());
    ASSERT(set.count(1) == 1, "big sphere intrudes into box");
    ASSERT(set.count(2) == 0, "small sphere does not reach box");
    ASSERT(set.count(3) == 1, "sphere fully inside box");
}

TEST(spatial_hash_mixed_points_and_spheres) {
    // Mixed inserts should not corrupt either query mode.
    SpatialHash3D h(1.0f);
    h.insert(Vec3{0, 0, 0}, 10);
    h.insert(Sphere{{3, 0, 0}, 2.5f}, 20); // reaches x=0.5
    h.insert(Vec3{5, 0, 0}, 30);

    std::vector<int32_t> hits;
    h.radiusQuery(Vec3{0, 0, 0}, 0.6f, hits);
    std::set<int32_t> set(hits.begin(), hits.end());
    ASSERT(set.count(10) == 1, "point at origin matches");
    ASSERT(set.count(20) == 1, "sphere reaching origin matches");
    ASSERT(set.count(30) == 0, "distant point does not match");

    hits.clear();
    AABB3 box{{4, -1, -1}, {6, 1, 1}};
    h.queryAABB(box, hits);
    std::set<int32_t> set2(hits.begin(), hits.end());
    ASSERT(set2.count(30) == 1, "point in AABB matches");
    ASSERT(set2.count(20) == 1, "sphere intrudes into AABB");
}
