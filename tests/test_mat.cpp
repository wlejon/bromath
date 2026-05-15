#include "test_framework.h"
#include "bromath/mat.h"
#include "bromath/transform.h"

using namespace bromath;

TEST(mat_identity) {
    Mat4 m = midentity();
    ASSERT(m.at(0,0) == 1.0f && m.at(1,1) == 1.0f && m.at(2,2) == 1.0f, "identity diag");
    ASSERT(m.at(0,1) == 0.0f && m.at(1,0) == 0.0f, "identity off-diag");
    Vec3 p = mtransformPoint(m, Vec3{3, 4, 5});
    ASSERT(p.x == 3 && p.y == 4 && p.z == 5, "identity * point");
}

TEST(mat_translate) {
    Mat4 m = mtranslate(Vec3{10, 20, 30});
    Vec3 p = mtransformPoint(m, Vec3{1, 2, 3});
    ASSERT(p.x == 11 && p.y == 22 && p.z == 33, "translate point");
    // direction ignores translation
    Vec3 d = mtransformDir(m, Vec3{1, 0, 0});
    ASSERT(d.x == 1 && d.y == 0 && d.z == 0, "translate dir unchanged");
}

TEST(mat_scale) {
    Mat4 m = mscale(Vec3{2, 3, 4});
    Vec3 p = mtransformPoint(m, Vec3{1, 1, 1});
    ASSERT(p.x == 2 && p.y == 3 && p.z == 4, "scale point");
}

TEST(mat_mul_associative_with_transform) {
    Mat4 t = mtranslate(Vec3{1, 2, 3});
    Mat4 s = mscale(Vec3{2, 2, 2});
    Mat4 ts = mmul(t, s);  // translate after scale
    Vec3 p = mtransformPoint(ts, Vec3{1, 1, 1});
    ASSERT(p.x == 3 && p.y == 4 && p.z == 5, "T*S point");
}

TEST(mat_transpose) {
    Mat4 m;
    m.at(0,1) = 5.0f;
    m.at(2,3) = 7.0f;
    Mat4 t = mtranspose(m);
    ASSERT(t.at(1,0) == 5.0f, "transpose 01->10");
    ASSERT(t.at(3,2) == 7.0f, "transpose 23->32");
}

TEST(mat_inverse) {
    Mat4 m = mmul(mtranslate(Vec3{1, 2, 3}), mscale(Vec3{2, 4, 8}));
    Mat4 inv = minverse(m);
    Mat4 id = mmul(m, inv);
    for (int i = 0; i < 16; ++i) {
        int row = i % 4, col = i / 4;
        float expected = (row == col) ? 1.0f : 0.0f;
        ASSERT(nearly(id.data[i], expected, 1e-4f), "m * m^-1 = I");
    }
}

TEST(mat_fromQuat) {
    Quat q = qaxisAngle(Vec3{0, 1, 0}, HALF_PI);
    Mat4 m = mfromQuat(q);
    Vec3 r = mtransformDir(m, Vec3{1, 0, 0});
    ASSERT(nearly(r.x, 0.0f, 1e-4f) && nearly(r.z, -1.0f, 1e-4f), "mat from Y-90");
}

TEST(mat_fromTRS_decompose_roundtrip) {
    Vec3 t{1, 2, 3};
    Quat r = qnorm(qaxisAngle(Vec3{0, 1, 0}, 0.7f));
    Vec3 s{2, 2, 2};
    Mat4 m = mfromTRS(t, r, s);
    Vec3 t2; Quat r2; Vec3 s2;
    mdecompose(m, t2, r2, s2);
    ASSERT(nearly(t2.x, t.x, 1e-4f) && nearly(t2.y, t.y, 1e-4f) && nearly(t2.z, t.z, 1e-4f), "T roundtrip");
    ASSERT(nearly(s2.x, s.x, 1e-4f) && nearly(s2.y, s.y, 1e-4f) && nearly(s2.z, s.z, 1e-4f), "S roundtrip");
    // Quaternion sign ambiguity — compare via rotated basis vector.
    Vec3 vr1 = qrotate(r, Vec3{1, 0, 0});
    Vec3 vr2 = qrotate(r2, Vec3{1, 0, 0});
    ASSERT(nearly(vr1.x, vr2.x, 1e-4f) && nearly(vr1.z, vr2.z, 1e-4f), "R roundtrip via basis");
}

TEST(mat_lookAt) {
    // Eye at +5z looking at origin, up = +y.
    Mat4 v = mlookAt(Vec3{0, 0, 5}, Vec3{0, 0, 0}, Vec3{0, 1, 0});
    // Origin in view space is in front of camera at -5z (right-handed).
    Vec3 p = mtransformPoint(v, Vec3{0, 0, 0});
    ASSERT(nearly(p.z, -5.0f, 1e-4f), "lookAt origin -> -5z");
}

TEST(mat_perspective_ortho) {
    Mat4 p = mperspective(deg2rad(90.0f), 1.0f, 0.1f, 100.0f);
    // 45 = fov/2 -> tan(45) = 1, so m00 = 1/aspect = 1.
    ASSERT(nearly(p.at(0,0), 1.0f, 1e-4f), "perspective 90deg, aspect 1");
    Mat4 o = mortho(-1, 1, -1, 1, 0.1f, 100.0f);
    ASSERT(nearly(o.at(0,0), 1.0f), "ortho x-scale");
}

TEST(transform_compose) {
    Transform parent;
    parent.position = {10, 0, 0};
    parent.rotation = qaxisAngle(Vec3{0, 1, 0}, HALF_PI);

    Transform child;
    child.position = {1, 0, 0};
    Transform combined = tmul(parent, child);
    // Parent's Y-90 turns child's +X offset into -Z, then add parent pos.
    ASSERT(nearly(combined.position.x, 10.0f, 1e-4f), "compose x");
    ASSERT(nearly(combined.position.z, -1.0f, 1e-4f), "compose z");
}

TEST(transform_mat4_roundtrip) {
    Transform t;
    t.position = {1, 2, 3};
    t.rotation = qnorm(qaxisAngle(Vec3{0, 1, 0}, 0.5f));
    t.scale    = {2, 2, 2};
    Mat4 m = ttoMat4(t);
    Transform back = tfromMat4(m);
    ASSERT(nearly(back.position.x, 1.0f, 1e-4f), "mat4 roundtrip pos");
    ASSERT(nearly(back.scale.x, 2.0f, 1e-4f), "mat4 roundtrip scale");
}
