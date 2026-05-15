#include "test_framework.h"
#include "bromath/scalar.h"
#include "bromath/angle.h"

using namespace bromath;

TEST(scalar_constants) {
    ASSERT(nearly(PI, 3.14159265f), "PI");
    ASSERT(nearly(TWO_PI, 2.0f * PI), "TWO_PI");
    ASSERT(nearly(HALF_PI, PI * 0.5f), "HALF_PI");
    ASSERT(nearly(deg2rad(180.0f), PI), "deg2rad");
    ASSERT(nearly(rad2deg(PI), 180.0f), "rad2deg");
}

TEST(scalar_clamp_saturate) {
    ASSERT(clamp(5.0f, 0.0f, 10.0f) == 5.0f, "clamp inside");
    ASSERT(clamp(-1.0f, 0.0f, 10.0f) == 0.0f, "clamp lo");
    ASSERT(clamp(11.0f, 0.0f, 10.0f) == 10.0f, "clamp hi");
    ASSERT(saturate(-0.5f) == 0.0f, "saturate lo");
    ASSERT(saturate(1.5f) == 1.0f, "saturate hi");
    ASSERT(saturate(0.3f) == 0.3f, "saturate mid");
}

TEST(scalar_lerp_remap) {
    ASSERT(nearly(lerp(0.0f, 10.0f, 0.5f), 5.0f), "lerp mid");
    ASSERT(nearly(mix(0.0f, 10.0f, 0.25f), 2.5f), "mix");
    ASSERT(nearly(invLerp(0.0f, 10.0f, 5.0f), 0.5f), "invLerp");
    ASSERT(nearly(remap(5.0f, 0.0f, 10.0f, 100.0f, 200.0f), 150.0f), "remap");
}

TEST(scalar_smoothstep) {
    ASSERT(smoothstep(0.0f, 1.0f, -1.0f) == 0.0f, "smoothstep below");
    ASSERT(smoothstep(0.0f, 1.0f, 2.0f) == 1.0f, "smoothstep above");
    ASSERT(nearly(smoothstep(0.0f, 1.0f, 0.5f), 0.5f), "smoothstep mid");
    ASSERT(nearly(smoothstep01(0.5f), 0.5f), "smoothstep01 mid");
    ASSERT(nearly(smootherstep(0.0f, 1.0f, 0.5f), 0.5f), "smootherstep mid");
}

TEST(scalar_misc) {
    ASSERT(sign(-3.0f) == -1.0f, "sign neg");
    ASSERT(sign(3.0f) == 1.0f, "sign pos");
    ASSERT(sign(0.0f) == 0.0f, "sign zero");
    ASSERT(abs(-2.5f) == 2.5f, "abs");
    ASSERT(sqr(4.0f) == 16.0f, "sqr");
    ASSERT(step(0.5f, 0.3f) == 0.0f, "step below");
    ASSERT(step(0.5f, 0.7f) == 1.0f, "step above");
    ASSERT(min(2, 5) == 2, "min");
    ASSERT(max(2, 5) == 5, "max");
}

TEST(angle_wrap) {
    ASSERT(nearly(wrapAngle(0.0f), 0.0f), "wrap 0");
    ASSERT(nearly(wrapAngle(PI), PI, 1e-4f) || nearly(wrapAngle(PI), -PI, 1e-4f), "wrap PI edge");
    ASSERT(nearly(wrapAngle(TWO_PI + 0.5f), 0.5f, 1e-4f), "wrap 2PI+x");
    ASSERT(nearly(wrapAngle(-TWO_PI - 0.5f), -0.5f, 1e-4f), "wrap -2PI-x");
    ASSERT(wrapAngle2Pi(-0.1f) > 0.0f, "wrap2Pi negative -> positive");
}

TEST(angle_delta_lerp) {
    ASSERT(nearly(angleDelta(0.0f, HALF_PI), HALF_PI), "delta 0->π/2");
    // -π to +π should choose short way (~0)
    float d = angleDelta(-PI + 0.1f, PI - 0.1f);
    ASSERT(std::fabs(d) < 0.3f, "delta wraps short way");
    ASSERT(nearly(angleLerp(0.0f, HALF_PI, 0.5f), HALF_PI * 0.5f), "angleLerp mid");
}
