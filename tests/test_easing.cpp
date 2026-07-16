#include "test_framework.h"
#include "bromath/easing.h"

using namespace bromath;

TEST(easing_endpoints) {
    // Every curve must satisfy f(0) = 0 and f(1) = 1.
    using Fn = float(*)(float);
    Fn fns[] = {
        easeLinear,
        easeQuadIn, easeQuadOut, easeQuadInOut,
        easeCubicIn, easeCubicOut, easeCubicInOut,
        easeQuartIn, easeQuartOut, easeQuartInOut,
        easeQuintIn, easeQuintOut, easeQuintInOut,
        easeSineIn, easeSineOut, easeSineInOut,
        easeExpoIn, easeExpoOut, easeExpoInOut,
        easeCircIn, easeCircOut, easeCircInOut,
        easeBackIn, easeBackOut, easeBackInOut,
        easeElasticIn, easeElasticOut, easeElasticInOut,
        easeBounceIn, easeBounceOut, easeBounceInOut,
    };
    for (Fn f : fns) {
        ASSERT(nearly(f(0.0f), 0.0f, 1e-5f), "ease(0) = 0");
        ASSERT(nearly(f(1.0f), 1.0f, 1e-5f), "ease(1) = 1");
    }
}

TEST(easing_values) {
    ASSERT(nearly(easeQuadIn(0.5f), 0.25f), "quadIn(0.5)");
    ASSERT(nearly(easeQuadOut(0.5f), 0.75f), "quadOut(0.5)");
    ASSERT(nearly(easeQuadInOut(0.25f), 0.125f), "quadInOut(0.25)");
    ASSERT(nearly(easeQuadInOut(0.75f), 0.875f), "quadInOut(0.75)");
    ASSERT(nearly(easeCubicIn(0.5f), 0.125f), "cubicIn(0.5)");
    ASSERT(nearly(easeSineInOut(0.5f), 0.5f), "sineInOut symmetric at midpoint");
    ASSERT(nearly(easeExpoIn(0.5f), std::pow(2.0f, -5.0f)), "expoIn(0.5)");
    // back overshoots below 0 early, elastic oscillates past 1 late.
    ASSERT(easeBackIn(0.3f) < 0.0f, "backIn dips negative");
    ASSERT(easeElasticOut(0.15f) > 1.0f, "elasticOut overshoots 1");
    // bounce stays within [0,1].
    for (int i = 0; i <= 20; i++) {
        float v = easeBounceOut(i / 20.0f);
        ASSERT(v >= -1e-5f && v <= 1.0f + 1e-5f, "bounceOut in range");
    }
    // in/out symmetry: easeXIn(t) == 1 - easeXOut(1-t).
    for (int i = 0; i <= 10; i++) {
        float t = i / 10.0f;
        ASSERT(nearly(easeCubicIn(t), 1.0f - easeCubicOut(1.0f - t), 1e-5f),
               "cubic in/out mirrored");
        ASSERT(nearly(easeBounceIn(t), 1.0f - easeBounceOut(1.0f - t), 1e-5f),
               "bounce in/out mirrored");
    }
}
