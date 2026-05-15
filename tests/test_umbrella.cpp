// Sanity check that the umbrella header pulls everything in cleanly.
#include "test_framework.h"
#include "bromath/bromath.h"

using namespace bromath;

TEST(umbrella_compiles) {
    Vec3 v = vnorm(Vec3{1, 1, 1});
    Quat q = qaxisAngle(v, HALF_PI);
    Mat4 m = mfromQuat(q);
    AABB3 a{{0,0,0}, {1,1,1}};
    Frustum f = ffromViewProj(mmul(mperspective(deg2rad(60.0f), 1.0f, 0.1f, 100.0f),
                                    mlookAt(Vec3{0,0,5}, Vec3{0,0,0}, Vec3{0,1,0})));
    Color c = cfromHex("#FF8800");
    uint64_t seed = 1;
    float x = randFloat01(seed);
    Smoother s; smootherSetTime(s, 5.0f, 44100.0f);
    GridFootprint2D g; g.cellSize = 1.0f; g.width = 4; g.depth = 4;
    ASSERT(m.at(3,3) == 1.0f, "mat identity element");
    ASSERT(a.max.x == 1.0f, "aabb max");
    ASSERT(f.planes[0].normal.x != 0.0f || f.planes[0].normal.y != 0.0f
            || f.planes[0].normal.z != 0.0f, "frustum plane non-zero");
    ASSERT(c.r > 0.5f, "hex orange has red");
    ASSERT(x >= 0.0f && x < 1.0f, "rng in range");
    ASSERT(s.coeff > 0.0f, "smoother coeff set");
    ASSERT(gridIndex2D(g, 1, 1) == 5, "grid index");
}
