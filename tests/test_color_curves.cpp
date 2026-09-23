#include "test_framework.h"
#include "bromath/color.h"
#include "bromath/curves.h"

#include <algorithm>
#include <cmath>

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

TEST(color_parse_css) {
    uint8_t r = 1, g = 2, b = 3, a = 4;
    // A literal must resolve (it is ambiguous between string/string_view
    // without the const char* overload).
    ASSERT(parseCSSColor("#fa0", r, g, b, a) && r == 255 && g == 170 && b == 0 && a == 255, "#RGB");
    ASSERT(parseCSSColor("#fa08", r, g, b, a) && r == 255 && g == 170 && b == 0 && a == 136, "#RGBA");
    ASSERT(parseCSSColor("#1A2b3C", r, g, b, a) && r == 0x1A && g == 0x2B && b == 0x3C && a == 255, "#RRGGBB");
    ASSERT(parseCSSColor(std::string("#11223380"), r, g, b, a) && r == 0x11 && b == 0x33 && a == 0x80, "#RRGGBBAA");
    ASSERT(parseCSSColor(std::string_view("Orange"), r, g, b, a) && r == 255 && g == 165 && b == 0, "names are case-insensitive");

    r = 7; g = 8; b = 9;
    ASSERT(!parseCSSColor("#GG0000", r, g, b, a), "non-hex digits are rejected, not read as 0");
    ASSERT(!parseCSSColor("#12345", r, g, b, a), "5 digits");
    ASSERT(!parseCSSColor("#", r, g, b, a), "bare #");
    ASSERT(!parseCSSColor("", r, g, b, a), "empty");
    ASSERT(!parseCSSColor("notacolor", r, g, b, a), "unknown name");
    ASSERT(!parseCSSColor(static_cast<const char*>(nullptr), r, g, b, a), "null");
    ASSERT(r == 7 && g == 8 && b == 9, "a failed parse leaves the channels alone");

    ASSERT(serializeCSSColor(0x1A, 0x2B, 0x3C) == "#1A2B3C", "serialize opaque");
    ASSERT(serializeCSSColor(Color8{1, 2, 3, 4}) == "#01020304", "serialize with alpha");
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

TEST(curve_catmullRom_repeated_ends) {
    // A spline starting on a repeated point stays a curve (it used to fall
    // back to a straight lerp for the whole segment).
    Vec3 p1{0, 0, 0}, p2{1, 0, 0}, p3{2, 1, 0};
    Vec3 mid = ccatmullRom(p1, p1, p2, p3, 0.5f);
    ASSERT(std::fabs(mid.y) > 1e-3f, "repeated first point: the segment still bends");
    ASSERT(nearly(ccatmullRom(p1, p1, p2, p3, 0.0f).x, 0.0f, 1e-4f) &&
               nearly(ccatmullRom(p1, p1, p2, p3, 1.0f).x, 1.0f, 1e-4f),
           "repeated first point: endpoints kept");
    Vec3 endMid = ccatmullRom(Vec3{-1, 1, 0}, p1, p2, p2, 0.5f);
    ASSERT(std::fabs(endMid.y) > 1e-3f, "repeated last point: the segment still bends");
    Vec3 same = ccatmullRom(Vec3{-1, 0, 0}, p1, p1, p3, 0.3f);
    ASSERT(same.x == 0 && same.y == 0 && same.z == 0, "zero-length segment stays at its point");
}

// y for x on the CSS cubic-bezier, by bisection in double (the reference).
static double easeRef(const CubicEase& c, double x) {
    auto bx = [&](double t) { double u = 1 - t; return 3*u*u*t*c.p1x + 3*u*t*t*c.p2x + t*t*t; };
    double lo = 0, hi = 1, t = 0.5;
    for (int i = 0; i < 100; ++i) {
        t = 0.5 * (lo + hi);
        if (bx(t) < x) lo = t; else hi = t;
    }
    double u = 1 - t;
    return 3*u*u*t*c.p1y + 3*u*t*t*c.p2y + t*t*t;
}

TEST(curve_cubicEase_flat_x) {
    // x(t) flat at t = 0.5 (control x = 1 then 0): Newton creeps toward the
    // inflection and stopped far short after 8 steps.
    CubicEase c{1.0f, 0.0f, 0.0f, 1.0f};
    float worst = 0.0f;
    for (float x : {0.3f, 0.45f, 0.49f, 0.499f, 0.501f, 0.51f, 0.55f, 0.7f}) {
        worst = std::max(worst, std::fabs(ccubicEase(c, x) - static_cast<float>(easeRef(c, x))));
    }
    ASSERT(worst < 2e-5f, "cubic ease matches the reference where x(t) is flat");
    CubicEase ease{0.25f, 0.1f, 0.25f, 1.0f};
    float worstEase = 0.0f;
    for (int i = 1; i < 100; ++i) {
        float x = i / 100.0f;
        worstEase = std::max(worstEase,
                             std::fabs(ccubicEase(ease, x) - static_cast<float>(easeRef(ease, x))));
    }
    ASSERT(worstEase < 1e-4f, "CSS ease matches the reference");
}

TEST(color_hex_trailing_and_nan) {
    Color bad = cfromHex("#FF000080Z");
    ASSERT(bad.a == 0.0f, "trailing characters after 8 digits are rejected");
    ASSERT(cfromHex("#FF000080").a > 0.49f, "exactly 8 digits still parse");
    Color8 n = ctoColor8(Color{std::nanf(""), 0.5f, 2.0f, std::nanf("")});
    ASSERT(n.r == 0 && n.b == 255 && n.a == 0, "NaN channels encode as 0");
}
