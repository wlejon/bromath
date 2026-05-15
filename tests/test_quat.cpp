#include "test_framework.h"
#include "bromath/quat.h"

using namespace bromath;

TEST(quat_identity) {
    Quat q = qidentity();
    ASSERT(q.x == 0 && q.y == 0 && q.z == 0 && q.w == 1, "identity");
    Vec3 v{1, 2, 3};
    Vec3 r = qrotate(q, v);
    ASSERT(nearly(r.x, 1.0f) && nearly(r.y, 2.0f) && nearly(r.z, 3.0f), "identity rotates to self");
}

TEST(quat_axis_angle) {
    // 90° about Y rotates +X to -Z
    Quat q = qaxisAngle(Vec3{0, 1, 0}, HALF_PI);
    Vec3 r = qrotate(q, Vec3{1, 0, 0});
    ASSERT(nearly(r.x, 0.0f, 1e-4f) && nearly(r.z, -1.0f, 1e-4f), "Y-90 rotates +X to -Z");
}

TEST(quat_mul) {
    // Two 90° Y rotations = 180° rotation: +X -> -X
    Quat q90 = qaxisAngle(Vec3{0, 1, 0}, HALF_PI);
    Quat q180 = qmul(q90, q90);
    Vec3 r = qrotate(q180, Vec3{1, 0, 0});
    ASSERT(nearly(r.x, -1.0f, 1e-4f) && nearly(r.z, 0.0f, 1e-4f), "Y-90 * Y-90 = Y-180");
}

TEST(quat_inverse_conjugate) {
    Quat q = qaxisAngle(Vec3{0, 1, 0}, 0.7f);
    Quat qi = qinverse(q);
    Quat id = qmul(q, qi);
    ASSERT(nearly(id.w, 1.0f, 1e-5f), "q * q^-1 = identity (w)");
    ASSERT(nearly(id.x, 0.0f, 1e-5f), "q * q^-1 = identity (x)");
}

TEST(quat_fromTo) {
    Quat q = qfromTo(Vec3{1, 0, 0}, Vec3{0, 1, 0});
    Vec3 r = qrotate(q, Vec3{1, 0, 0});
    ASSERT(nearly(r.y, 1.0f, 1e-4f), "fromTo +X -> +Y");
    // 180 flip case
    Quat q180 = qfromTo(Vec3{1, 0, 0}, Vec3{-1, 0, 0});
    Vec3 r180 = qrotate(q180, Vec3{1, 0, 0});
    ASSERT(nearly(r180.x, -1.0f, 1e-4f), "fromTo +X -> -X (180)");
}

TEST(quat_euler_roundtrip) {
    Vec3 e{0.3f, 0.5f, -0.7f};
    Quat q = qfromEuler(e);
    Vec3 back = qtoEuler(qnorm(q));
    ASSERT(nearly(back.x, e.x, 1e-3f), "euler roundtrip x");
    ASSERT(nearly(back.y, e.y, 1e-3f), "euler roundtrip y");
    ASSERT(nearly(back.z, e.z, 1e-3f), "euler roundtrip z");
}

TEST(quat_slerp) {
    Quat a = qidentity();
    Quat b = qaxisAngle(Vec3{0, 1, 0}, HALF_PI);
    Quat half = qslerp(a, b, 0.5f);
    Vec3 r = qrotate(half, Vec3{1, 0, 0});
    // halfway between +X and -Z
    float expected = std::cos(HALF_PI * 0.5f);
    ASSERT(nearly(r.x, expected, 1e-3f), "slerp halfway x");
    ASSERT(nearly(r.z, -expected, 1e-3f), "slerp halfway z");
}

TEST(quat_norm) {
    Quat q{2, 0, 0, 0};
    Quat n = qnorm(q);
    ASSERT(nearly(qlen(n), 1.0f), "normalize length 1");
    ASSERT(nearly(n.x, 1.0f), "normalize x");
}
