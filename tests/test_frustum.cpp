#include "test_framework.h"
#include "bromath/frustum.h"

using namespace bromath;

static Mat4 makeVP() {
    // Camera at origin looking down -Z, 90° fov, 1:1 aspect, near=0.1 far=100.
    Mat4 view = mlookAt(Vec3{0, 0, 0}, Vec3{0, 0, -1}, Vec3{0, 1, 0});
    Mat4 proj = mperspective(deg2rad(90.0f), 1.0f, 0.1f, 100.0f);
    return mmul(proj, view);
}

TEST(frustum_contains_point) {
    Frustum f = ffromViewProj(makeVP());
    ASSERT(fcontains(f, Vec3{0, 0, -10}), "point in front of camera");
    ASSERT(!fcontains(f, Vec3{0, 0, 10}), "point behind camera");
    ASSERT(!fcontains(f, Vec3{0, 0, -200}), "point past far plane");
    ASSERT(!fcontains(f, Vec3{0, 0, -0.05f}), "point closer than near");
}

TEST(frustum_aabb_cull) {
    Frustum f = ffromViewProj(makeVP());
    AABB3 inside{{-1, -1, -11}, {1, 1, -9}};
    ASSERT(fintersects(f, inside), "AABB centered in front");
    AABB3 behind{{-1, -1, 10}, {1, 1, 12}};
    ASSERT(!fintersects(f, behind), "AABB fully behind camera");
    AABB3 far{{-1, -1, -200}, {1, 1, -150}};
    ASSERT(!fintersects(f, far), "AABB past far plane");
    AABB3 spanning{{-1, -1, -50}, {1, 1, 50}};
    ASSERT(fintersects(f, spanning), "AABB spanning frustum kept");
}

TEST(frustum_sphere_cull) {
    Frustum f = ffromViewProj(makeVP());
    Sphere visible{{0, 0, -10}, 1.0f};
    ASSERT(fintersects(f, visible), "sphere in front");
    Sphere behind{{0, 0, 10}, 1.0f};
    ASSERT(!fintersects(f, behind), "sphere behind camera");
    // Sphere just past far plane should be culled.
    Sphere tooFar{{0, 0, -200}, 1.0f};
    ASSERT(!fintersects(f, tooFar), "sphere past far");
    // Sphere overlapping far plane should be kept.
    Sphere onFar{{0, 0, -99}, 5.0f};
    ASSERT(fintersects(f, onFar), "sphere crossing far plane");
}
