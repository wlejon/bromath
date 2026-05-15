#include "test_framework.h"
#include "bromath/color.h"
#include "bromath/curves.h"

using namespace bromath;

TEST(color_lerp) {
    Color a{1, 0, 0, 1};
    Color b{0, 1, 0, 1};
    Color m = clerp(a, b, 0.5f);
    ASSERT(nearly(m.r, 0.5f) && nearly(m.g, 0.5f), "color lerp");
}

TEST(color_srgb_roundtrip) {
    // Mid-grey in sRGB should NOT be 0.5 in linear (gamma effect).
    Color8 grey{128, 128, 128, 255};
    Color lin = cfromColor8(grey);
    ASSERT(lin.r < 0.25f, "sRGB 128 -> linear < 0.25");
    Color8 back = ctoColor8(lin);
    ASSERT(back.r >= 127 && back.r <= 129, "roundtrip preserves byte");
}

TEST(color_from_hex) {
    Color red = cfromHex("#FF0000");
    ASSERT(nearly(red.r, 1.0f, 1e-3f), "hex red is full red");
    ASSERT(nearly(red.g, 0.0f), "hex red has no green");
    ASSERT(nearly(red.a, 1.0f), "hex red opaque");
    Color half = cfromHex("#0000FF80");
    ASSERT(nearly(half.a, 128.0f / 255.0f, 1e-3f), "hex alpha");
    Color bad = cfromHex("nope");
    ASSERT(bad.a == 0.0f, "bad hex returns transparent");
}

TEST(color_hsv) {
    // Pure red in HSV.
    Color red = cfromHSV(0.0f, 1.0f, 1.0f);
    ASSERT(nearly(red.r, 1.0f, 1e-3f), "HSV red");
    Color green = cfromHSV(120.0f, 1.0f, 1.0f);
    ASSERT(nearly(green.g, 1.0f, 1e-3f), "HSV green");
    Color blue = cfromHSV(240.0f, 1.0f, 1.0f);
    ASSERT(nearly(blue.b, 1.0f, 1e-3f), "HSV blue");
}

TEST(curve_cubicEase) {
    CubicEase ease{0.25f, 0.1f, 0.25f, 1.0f}; // CSS "ease"
    ASSERT(ccubicEase(ease, 0.0f) == 0.0f, "ease(0) = 0");
    ASSERT(ccubicEase(ease, 1.0f) == 1.0f, "ease(1) = 1");
    float mid = ccubicEase(ease, 0.5f);
    ASSERT(mid > 0.5f, "CSS ease accelerates early");
}

TEST(curve_bezier) {
    Vec3 p0{0,0,0}, p1{0,1,0}, p2{1,1,0}, p3{1,0,0};
    Vec3 s = cbezier(p0, p1, p2, p3, 0.0f);
    ASSERT(s.x == 0 && s.y == 0, "bezier(0) = p0");
    Vec3 e = cbezier(p0, p1, p2, p3, 1.0f);
    ASSERT(e.x == 1 && e.y == 0, "bezier(1) = p3");
    Vec3 m = cbezier(p0, p1, p2, p3, 0.5f);
    ASSERT(m.y > 0.5f, "midpoint pulled toward middle controls");
    Vec3 tan = cbezierTangent(p0, p1, p2, p3, 0.0f);
    ASSERT(tan.y > 0.0f, "start tangent points up toward p1");
}

TEST(curve_hermite) {
    Vec3 p0{0, 0, 0}, p1{10, 0, 0};
    Vec3 m0{0, 5, 0}, m1{0, -5, 0};
    Vec3 s = chermite(p0, m0, p1, m1, 0.0f);
    ASSERT(s.x == 0, "hermite(0) = p0");
    Vec3 e = chermite(p0, m0, p1, m1, 1.0f);
    ASSERT(nearly(e.x, 10.0f), "hermite(1) = p1");
    Vec3 mid = chermite(p0, m0, p1, m1, 0.5f);
    ASSERT(mid.y > 0.0f, "hermite arcs above straight line");
}

TEST(curve_catmullRom) {
    Vec3 p0{-1, 0, 0}, p1{0, 0, 0}, p2{1, 0, 0}, p3{2, 0, 0};
    Vec3 s = ccatmullRom(p0, p1, p2, p3, 0.0f);
    ASSERT(nearly(s.x, p1.x, 1e-3f), "catmull(0) = p1");
    Vec3 e = ccatmullRom(p0, p1, p2, p3, 1.0f);
    ASSERT(nearly(e.x, p2.x, 1e-3f), "catmull(1) = p2");
}
