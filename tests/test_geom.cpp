#include "test_framework.h"
#include "bromath/aabb.h"
#include "bromath/plane.h"
#include "bromath/sphere.h"
#include "bromath/ray.h"

using namespace bromath;

TEST(aabb_contains_intersect) {
    AABB3 a{{0, 0, 0}, {1, 1, 1}};
    ASSERT(acontains(a, Vec3{0.5f, 0.5f, 0.5f}), "contains inside");
    ASSERT(!acontains(a, Vec3{2, 0, 0}), "contains outside");
    AABB3 b{{0.5f, 0.5f, 0.5f}, {2, 2, 2}};
    ASSERT(aintersects(a, b), "intersects overlap");
    AABB3 c{{2, 2, 2}, {3, 3, 3}};
    ASSERT(!aintersects(a, c), "intersects disjoint");
}

TEST(aabb_expand_merge) {
    AABB3 a = aempty3();
    ASSERT(aisEmpty(a), "empty box is empty");
    a = aexpand(a, Vec3{1, 2, 3});
    a = aexpand(a, Vec3{-1, 5, 0});
    ASSERT(a.min.x == -1 && a.max.y == 5, "expand grows");
    AABB3 b{{10, 10, 10}, {20, 20, 20}};
    AABB3 m = amerge(a, b);
    ASSERT(m.max.x == 20 && m.min.x == -1, "merge spans both");
}

TEST(aabb_fromPoints) {
    float pts[] = {0,0,0, 1,2,3, -1,-2,-3, 5,0,0};
    AABB3 a = afromPoints(pts, 4);
    ASSERT(a.min.x == -1 && a.max.x == 5, "fromPoints x");
    ASSERT(a.min.y == -2 && a.max.y == 2, "fromPoints y");
}

TEST(aabb_transform) {
    AABB3 a{{0, 0, 0}, {1, 1, 1}};
    Mat4 m = mtranslate(Vec3{10, 20, 30});
    AABB3 t = atransform(a, m);
    ASSERT(t.min.x == 10 && t.max.x == 11, "translated AABB");
}

TEST(plane_distance_project) {
    Plane p = pfromPointNormal(Vec3{0, 5, 0}, Vec3{0, 1, 0});
    ASSERT(nearly(psignedDistance(p, Vec3{3, 7, -1}), 2.0f), "signed dist above");
    ASSERT(nearly(psignedDistance(p, Vec3{0, 4, 0}), -1.0f), "signed dist below");
    Vec3 proj = pproject(p, Vec3{3, 7, -1});
    ASSERT(nearly(proj.y, 5.0f), "projection lands on plane");
}

TEST(plane_fromPoints) {
    Plane p = pfromPoints(Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 0, 1});
    // Normal should be +Y or -Y depending on winding (CCW = +Y here).
    ASSERT(nearly(std::fabs(p.normal.y), 1.0f, 1e-4f), "plane normal Y-ish");
}

TEST(sphere_contains_intersect) {
    Sphere s{{0, 0, 0}, 2.0f};
    ASSERT(scontains(s, Vec3{1, 0, 0}), "contains inside");
    ASSERT(!scontains(s, Vec3{3, 0, 0}), "contains outside");
    Sphere t{{3, 0, 0}, 1.5f};
    ASSERT(sintersects(s, t), "spheres just touch");
    Sphere u{{10, 0, 0}, 1.0f};
    ASSERT(!sintersects(s, u), "spheres disjoint");
}

TEST(sphere_intersect_volume) {
    Sphere s{{0, 0, 0}, 1.0f};
    Sphere t{{10, 0, 0}, 1.0f};
    ASSERT(sintersectVolume(s, t) == 0.0f, "disjoint volume 0");
    Sphere u{{0.5f, 0, 0}, 1.0f};
    float v = sintersectVolume(s, u);
    ASSERT(v > 0.0f, "overlap volume positive");
    Sphere big{{0, 0, 0}, 5.0f};
    Sphere small{{0, 0, 0}, 1.0f};
    float fully = sintersectVolume(big, small);
    float expectedSmall = (4.0f / 3.0f) * PI;
    ASSERT(nearly(fully, expectedSmall, 1e-4f), "containment = small volume");
}

TEST(ray_aabb_hit) {
    Ray r{ Vec3{-5, 0.5f, 0.5f}, Vec3{1, 0, 0} };
    AABB3 a{{0, 0, 0}, {1, 1, 1}};
    RayHit h = rIntersectAABB(r, a);
    ASSERT(h.hit, "ray hits box");
    ASSERT(nearly(h.t, 5.0f, 1e-4f), "ray-aabb t = 5");
    ASSERT(nearly(h.normal.x, -1.0f), "hit normal -X");
}

TEST(ray_aabb_miss) {
    Ray r{ Vec3{-5, 5, 5}, Vec3{1, 0, 0} };
    AABB3 a{{0, 0, 0}, {1, 1, 1}};
    RayHit h = rIntersectAABB(r, a);
    ASSERT(!h.hit, "ray misses box");
}

TEST(ray_sphere) {
    Ray r{ Vec3{0, 0, -5}, Vec3{0, 0, 1} };
    Sphere s{{0, 0, 0}, 1.0f};
    RayHit h = rIntersectSphere(r, s);
    ASSERT(h.hit, "ray hits sphere");
    ASSERT(nearly(h.t, 4.0f, 1e-4f), "t = 4 (5 - radius)");
    ASSERT(nearly(h.normal.z, -1.0f, 1e-4f), "normal points back at ray");
}

TEST(ray_plane) {
    Ray r{ Vec3{0, 10, 0}, Vec3{0, -1, 0} };
    Plane p = pfromPointNormal(Vec3{0, 0, 0}, Vec3{0, 1, 0});
    RayHit h = rIntersectPlane(r, p);
    ASSERT(h.hit, "ray hits plane");
    ASSERT(nearly(h.t, 10.0f), "t = 10");
}

TEST(ray_triangle) {
    Ray r{ Vec3{0.25f, 0.25f, -1}, Vec3{0, 0, 1} };
    Vec3 v0{0, 0, 0}, v1{1, 0, 0}, v2{0, 1, 0};
    RayHit h = rIntersectTriangle(r, v0, v1, v2);
    ASSERT(h.hit, "ray hits triangle");
    ASSERT(nearly(h.t, 1.0f, 1e-4f), "t = 1");
    Ray miss{ Vec3{5, 5, -1}, Vec3{0, 0, 1} };
    RayHit hm = rIntersectTriangle(miss, v0, v1, v2);
    ASSERT(!hm.hit, "ray misses triangle");
}
