#include "test_framework.h"
#include "bromath/vec.h"

using namespace bromath;

TEST(vec3_basic) {
    Vec3 a{1, 2, 3};
    Vec3 b{4, 5, 6};
    Vec3 s = a + b;
    ASSERT(s.x == 5 && s.y == 7 && s.z == 9, "add");
    Vec3 d = b - a;
    ASSERT(d.x == 3 && d.y == 3 && d.z == 3, "sub");
    Vec3 m = a * 2.0f;
    ASSERT(m.x == 2 && m.y == 4 && m.z == 6, "scale");
    Vec3 n = -a;
    ASSERT(n.x == -1 && n.y == -2 && n.z == -3, "neg");
}

TEST(vec3_dot_cross) {
    ASSERT(nearly(vdot(Vec3{1,2,3}, Vec3{4,5,6}), 32.0f), "dot");
    Vec3 c = vcross(Vec3{1,0,0}, Vec3{0,1,0});
    ASSERT(nearly(c.x, 0.0f) && nearly(c.y, 0.0f) && nearly(c.z, 1.0f), "cross x*y=z");
}

TEST(vec3_length_norm) {
    ASSERT(nearly(vlen(Vec3{3, 4, 0}), 5.0f), "length 3-4-5");
    ASSERT(nearly(vlen2(Vec3{1, 2, 2}), 9.0f), "length2");
    Vec3 n = vnorm(Vec3{0, 0, 5});
    ASSERT(nearly(n.z, 1.0f) && nearly(n.x, 0.0f), "normalize");
    Vec3 z = vnorm(Vec3{0, 0, 0});
    ASSERT(z.x == 0 && z.y == 0 && z.z == 0, "normalize zero -> zero");
    Vec3 f = vnormOr(Vec3{0,0,0}, Vec3{1,0,0});
    ASSERT(f.x == 1.0f, "normalizeOr fallback");
}

TEST(vec3_distance) {
    ASSERT(nearly(vdist(Vec3{0,0,0}, Vec3{3,4,0}), 5.0f), "dist");
    ASSERT(nearly(vdist2(Vec3{0,0,0}, Vec3{3,4,0}), 25.0f), "dist2");
}

TEST(vec3_lerp_reflect_project) {
    Vec3 m = vlerp(Vec3{0,0,0}, Vec3{10,20,30}, 0.5f);
    ASSERT(m.x == 5 && m.y == 10 && m.z == 15, "lerp");
    // reflect (1,-1,0) about up normal (0,1,0) -> (1,1,0)
    Vec3 r = vreflect(Vec3{1, -1, 0}, Vec3{0, 1, 0});
    ASSERT(nearly(r.x, 1.0f) && nearly(r.y, 1.0f), "reflect");
    Vec3 p = vproject(Vec3{3, 4, 0}, Vec3{1, 0, 0});
    ASSERT(nearly(p.x, 3.0f) && nearly(p.y, 0.0f), "project on x");
}

TEST(vec3_min_max_perp) {
    Vec3 lo = vmin(Vec3{1, 5, -3}, Vec3{2, 4, 0});
    ASSERT(lo.x == 1 && lo.y == 4 && lo.z == -3, "vmin");
    Vec3 hi = vmax(Vec3{1, 5, -3}, Vec3{2, 4, 0});
    ASSERT(hi.x == 2 && hi.y == 5 && hi.z == 0, "vmax");
    Vec3 perp = vperpendicular(Vec3{1, 0, 0});
    ASSERT(nearly(vdot(perp, Vec3{1,0,0}), 0.0f, 1e-5f), "perpendicular orthogonal");
}

TEST(vec2_basic) {
    Vec2 a{3, 4};
    ASSERT(nearly(vlen(a), 5.0f), "vec2 len");
    Vec2 n = vnorm(a);
    ASSERT(nearly(vlen(n), 1.0f), "vec2 normalized");
    // 2D cross: x cross y = +1
    ASSERT(nearly(vcross(Vec2{1,0}, Vec2{0,1}), 1.0f), "vec2 cross");
}
