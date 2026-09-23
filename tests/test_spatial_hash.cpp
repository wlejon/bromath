#include "test_framework.h"

#include "bromath/aabb.h"
#include "bromath/spatial_hash.h"
#include "bromath/sphere.h"
#include "bromath/vec.h"

#include <algorithm>
#include <cmath>
#include <limits>
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

TEST(spatial_hash_stress_insert_remove_many) {
    SpatialHash3D h(0.5f);
    const int N = 2000;
    std::mt19937 rng(12345);
    std::uniform_real_distribution<float> dist(-20.0f, 20.0f);
    std::uniform_real_distribution<float> rdist(0.1f, 1.5f);

    struct Item {
        Vec3 p;
        float r;
        int32_t id;
    };
    std::vector<Item> items;
    items.reserve(N);

    for (int i = 0; i < N; ++i) {
        Vec3 p{dist(rng), dist(rng), dist(rng)};
        float r = (i % 2 == 0) ? 0.0f : rdist(rng);
        int32_t id = i + 100;
        items.push_back({p, r, id});
        if (r > 0.0f) {
            h.insert(Sphere{p, r}, id);
        } else {
            h.insert(p, id);
        }
    }
    ASSERT(h.size() == static_cast<size_t>(N), "all items inserted");

    // Remove half the items
    for (int i = 0; i < N; i += 2) {
        h.remove(items[i].id);
    }
    ASSERT(h.size() == static_cast<size_t>(N / 2), "half items removed");

    // Verify remaining items are found accurately
    Vec3 qCenter{0, 0, 0};
    float qRadius = 10.0f;
    std::set<int32_t> brute;
    for (int i = 1; i < N; i += 2) {
        float reach = qRadius + items[i].r;
        if (vdist(items[i].p, qCenter) <= reach) {
            brute.insert(items[i].id);
        }
    }
    std::vector<int32_t> got;
    h.radiusQuery(qCenter, qRadius, got);
    std::set<int32_t> gotSet(got.begin(), got.end());
    ASSERT(brute == gotSet, "post-removal radius query matches ground truth");
}

TEST(spatial_hash_negative_and_zero_radius) {
    SpatialHash3D h(1.0f);
    h.insert(Vec3{0, 0, 0}, 1);
    h.insert(Sphere{{0.5f, 0, 0}, 3.0f}, 2);
    std::vector<int32_t> hits;
    h.radiusQuery(Vec3{0, 0, 0}, -5.0f, hits);
    ASSERT(hits.empty(), "negative radius matches nothing");
    h.radiusQuery(Vec3{0, 0, 0}, std::nanf(""), hits);
    ASSERT(hits.empty(), "NaN radius matches nothing");

    SpatialHash3D p(1.0f);
    p.insert(Vec3{2, 2, 2}, 7);
    p.radiusQuery(Vec3{2, 2, 2}, 0.0f, hits);
    ASSERT(hits.size() == 1 && hits[0] == 7, "radius 0 finds a point at the center");
}

TEST(spatial_hash_huge_query_no_duplicates) {
    // A query box wider than the 2^21-cell key period used to visit wrapped
    // cells more than once (reporting ids twice), and a huge radius looped
    // over every cell in its cube.
    SpatialHash3D h(1.0f);
    h.insert(Vec3{0, 0, 0}, 1);
    h.insert(Vec3{3e6f, 0, 0}, 2);
    h.insert(Vec3{-3e6f, 5, 5}, 3);
    std::vector<int32_t> hits;
    h.radiusQuery(Vec3{0, 0, 0}, 1e7f, hits);
    std::vector<int32_t> sorted = hits;
    std::sort(sorted.begin(), sorted.end());
    ASSERT(sorted == (std::vector<int32_t>{1, 2, 3}), "huge radius finds each entry once");

    hits.clear();
    h.queryAABB(AABB3{{-1e7f, -1e7f, -1e7f}, {1e7f, 1e7f, 1e7f}}, hits);
    ASSERT(hits.size() == 3, "huge AABB reports each entry once");

    ASSERT(h.nearest(Vec3{2.9e6f, 0, 0}, 1e9f) == 2, "nearest with a huge radius");

    // Infinite radius: everything, once.
    hits.clear();
    h.radiusQuery(Vec3{0, 0, 0}, std::numeric_limits<float>::infinity(), hits);
    ASSERT(hits.size() == 3, "infinite radius reports each entry once");
}

TEST(spatial_hash_nonfinite_positions) {
    // Out-of-int-range and NaN positions have a cell rather than an
    // undefined conversion; they never match a finite query.
    SpatialHash3D h(0.001f);
    h.insert(Vec3{1e30f, 0, 0}, 1);
    h.insert(Vec3{std::nanf(""), 0, 0}, 2);
    h.insert(Vec3{-std::numeric_limits<float>::infinity(), 0, 0}, 3);
    h.insert(Vec3{0, 0, 0}, 4);
    std::vector<int32_t> hits;
    h.radiusQuery(Vec3{0, 0, 0}, 0.01f, hits);
    ASSERT(hits.size() == 1 && hits[0] == 4, "only the finite nearby point matches");
    h.remove(2);
    h.remove(1);
    ASSERT(h.size() == 2, "non-finite entries remove cleanly");
}

TEST(spatial_hash_cell_walk_matches_scan) {
    // Small queries over a dense grid walk cells; large ones scan. Both must
    // agree with brute force.
    SpatialHash3D h(0.25f);
    std::mt19937 rng(7);
    std::uniform_real_distribution<float> uni(-5.0f, 5.0f);
    std::vector<Vec3> pts;
    for (int i = 0; i < 5000; ++i) {
        Vec3 p{uni(rng), uni(rng), uni(rng)};
        pts.push_back(p);
        h.insert(p, i);
    }
    for (float r : {0.1f, 0.6f, 2.0f, 20.0f}) {
        Vec3 c{0.3f, -0.2f, 0.1f};
        std::vector<int32_t> got;
        h.radiusQuery(c, r, got);
        std::vector<int32_t> brute;
        for (size_t i = 0; i < pts.size(); ++i)
            if (vdist2(pts[i], c) <= r * r) brute.push_back(static_cast<int32_t>(i));
        std::sort(got.begin(), got.end());
        ASSERT(got == brute, "radius query equals brute force at every scale");

        AABB3 box{{c.x - r, c.y - r, c.z - r}, {c.x + r, c.y + r, c.z + r}};
        std::vector<int32_t> gotBox;
        h.queryAABB(box, gotBox);
        std::vector<int32_t> bruteBox;
        for (size_t i = 0; i < pts.size(); ++i)
            if (acontains(box, pts[i])) bruteBox.push_back(static_cast<int32_t>(i));
        std::sort(gotBox.begin(), gotBox.end());
        ASSERT(gotBox == bruteBox, "AABB query equals brute force at every scale");
    }
}

