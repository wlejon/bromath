#include "test_framework.h"
#include "bromath/mat.h"
#include "bromath/transform.h"

#include <cmath>

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

TEST(mat_inverse_scale_relative) {
    // A well-conditioned matrix at a tiny scale: det = 1e-24, far under the
    // old fixed 1e-20 cutoff, but the matrix is invertible.
    const float scales[] = {1e-6f, 1e-3f, 1.0f, 1e3f, 1e6f};
    for (float k : scales) {
        Mat4 m = mfromTRS(Vec3{1, -2, 3},
                          qnorm(qaxisAngle(Vec3{1, 1, 0}, 0.6f)),
                          Vec3{k, 2 * k, 0.5f * k});
        Mat4 id = mmul(m, minverse(m));
        for (int i = 0; i < 16; ++i) {
            int row = i % 4, col = i / 4;
            float expected = (row == col) ? 1.0f : 0.0f;
            ASSERT(nearly(id.data[i], expected, 1e-3f), "m * m^-1 = I at any scale");
        }
    }
    // Anisotropic: one axis scaled 1e-5 against a large translation.
    Mat4 thin = mmul(mtranslate(Vec3{1e4f, 0, 0}), mscale(Vec3{1e-5f, 1, 1}));
    Mat4 thinInv = minverse(thin);
    ASSERT(nearly(thinInv.at(0, 0), 1e5f, 1.0f), "thin axis inverts");

    // Singular at any scale: a zero column, a repeated column, rank 2 scaled
    // tiny, and a NaN entry all give the identity.
    auto isIdentity = [](const Mat4& m) {
        for (int i = 0; i < 16; ++i)
            if (m.data[i] != ((i % 5 == 0) ? 1.0f : 0.0f)) return false;
        return true;
    };
    ASSERT(isIdentity(minverse(mscale(Vec3{1, 0, 1}))), "zero column");
    for (float k : scales) {
        Mat4 dup;
        for (int r = 0; r < 4; ++r) {
            dup.at(r, 0) = (r + 1) * 0.37f * k;
            dup.at(r, 1) = (r + 1) * 0.37f * k;   // == column 0
            dup.at(r, 2) = (3 - r) * 0.11f * k;
            dup.at(r, 3) = (r * r + 1) * 0.5f * k;
        }
        ASSERT(isIdentity(minverse(dup)), "repeated column is singular at any scale");
        Mat4 rank;  // rows 1,2,3 arithmetic progression: 4x4 of i+j form
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c) rank.at(r, c) = float(r * 4 + c + 1) * 0.1f * k;
        ASSERT(isIdentity(minverse(rank)), "rank-2 matrix is singular at any scale");
    }
    Mat4 bad = midentity();
    bad.at(1, 2) = NAN;
    ASSERT(isIdentity(minverse(bad)), "NaN entry");
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

static bool matNear(const Mat4& a, const Mat4& b, float eps = 1e-4f) {
    for (int i = 0; i < 16; ++i)
        if (!nearly(a.data[i], b.data[i], eps)) return false;
    return true;
}

TEST(mat_decompose_mirror) {
    // A mirrored TRS: the quaternion used to come out of an improper basis,
    // so rebuilding the matrix gave something else.
    Quat q = qnorm(Quat{0.2f, 0.4f, -0.1f, 0.9f});
    for (Vec3 s : {Vec3{-2, 1, 3}, Vec3{2, -1, 3}, Vec3{2, 1, -3}, Vec3{-1, -1, -1}}) {
        Mat4 m = mfromTRS(Vec3{1, 2, 3}, q, s);
        Vec3 t, sc;
        Quat r;
        mdecompose(m, t, r, sc);
        ASSERT(nearly(qlen(r), 1.0f, 1e-4f), "decomposed rotation is unit length");
        ASSERT(matNear(mfromTRS(t, r, sc), m), "mirrored TRS round-trips through decompose");
    }
    Vec3 t, sc;
    Quat r;
    mdecompose(mfromTRS(Vec3{}, q, Vec3{1, 2, 3}), t, r, sc);
    ASSERT(sc.x > 0 && sc.y > 0 && sc.z > 0, "a proper matrix keeps positive scales");
}

TEST(mat_look_at_parallel_up) {
    // Looking straight down with up = +y: right = forward x up is zero.
    Mat4 v = mlookAt(Vec3{0, 10, 0}, Vec3{0, 0, 0}, Vec3{0, 1, 0});
    bool finite = true;
    for (float x : v.data) finite = finite && std::isfinite(x);
    ASSERT(finite, "look-at along up is finite");
    Vec3 p = mtransformPoint(v, Vec3{0, 0, 0});
    ASSERT(nearly(p.x, 0.0f, 1e-4f) && nearly(p.y, 0.0f, 1e-4f) && nearly(p.z, -10.0f, 1e-4f),
           "look-at along up still puts the target 10 units ahead");
    // The basis stays orthonormal (non-degenerate).
    Vec3 r0{v.at(0, 0), v.at(0, 1), v.at(0, 2)};
    Vec3 r1{v.at(1, 0), v.at(1, 1), v.at(1, 2)};
    ASSERT(nearly(vlen(r0), 1.0f, 1e-4f) && nearly(vlen(r1), 1.0f, 1e-4f) &&
               nearly(vdot(r0, r1), 0.0f, 1e-4f),
           "look-at along up has an orthonormal basis");
    Mat4 z = mlookAt(Vec3{0, 0, 0}, Vec3{0, 0, -1}, Vec3{0, 0, 0});
    ASSERT(nearly(vlen(Vec3{z.at(0, 0), z.at(0, 1), z.at(0, 2)}), 1.0f, 1e-4f),
           "look-at with a zero up still has a right vector");
}
