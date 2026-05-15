#include "test_framework.h"
#include "bromath/smoother.h"
#include "bromath/grid.h"

using namespace bromath;

TEST(smoother_chase) {
    Smoother s;
    smootherSetTime(s, 10.0f, 44100.0f);
    smootherReset(s, 0.0f);
    smootherTarget(s, 1.0f);
    // After ~10ms (441 samples) should be ~95% of the way.
    float v = smootherTickN(s, 441);
    ASSERT(v > 0.9f && v < 1.0f, "smoother reaches ~95% in stated time");
    // After many more ticks, should converge.
    v = smootherTickN(s, 10000);
    ASSERT(nearly(v, 1.0f, 1e-4f), "smoother converges");
}

TEST(smoother_reset_zero_time) {
    Smoother s;
    smootherSetTime(s, 0.0f, 44100.0f);
    smootherReset(s, 0.0f);
    smootherTarget(s, 1.0f);
    float v = smootherTick(s);
    ASSERT(v == 1.0f, "zero time snaps to target");
}

TEST(grid_bounds_and_index) {
    GridFootprint2D g;
    g.origin = {0, 0};
    g.cellSize = 1.0f;
    g.width = 4;
    g.depth = 3;
    ASSERT(gridIndex2D(g, 2, 1) == 6, "row-major index");
    ASSERT(gridInBounds(g, 0, 0), "in bounds origin");
    ASSERT(!gridInBounds(g, 4, 0), "out of bounds x");
    ASSERT(!gridInBounds(g, 0, 3), "out of bounds y");
}

TEST(grid_cell_of) {
    GridFootprint2D g;
    g.origin = {10, 20};
    g.cellSize = 2.0f;
    g.width = 5;
    g.depth = 5;
    int c, r;
    gridCellOf(g, Vec2{10.5f, 20.5f}, c, r);
    ASSERT(c == 0 && r == 0, "cell origin");
    gridCellOf(g, Vec2{14.0f, 25.0f}, c, r);
    ASSERT(c == 2 && r == 2, "cell offset");
    Vec2 ctr = gridCellCenter(g, 0, 0);
    ASSERT(nearly(ctr.x, 11.0f) && nearly(ctr.y, 21.0f), "cell center");
}

TEST(grid_index_3d) {
    // 2x3x4 grid (width=2, height=3, depth=4) — x outer, y middle, z inner.
    int idx000 = gridIndex3D(0, 0, 0, 4, 3);
    int idx001 = gridIndex3D(0, 0, 1, 4, 3);
    int idx010 = gridIndex3D(0, 1, 0, 4, 3);
    int idx100 = gridIndex3D(1, 0, 0, 4, 3);
    ASSERT(idx000 == 0, "3D origin");
    ASSERT(idx001 == 1, "3D z stride 1");
    ASSERT(idx010 == 4, "3D y stride = depth");
    ASSERT(idx100 == 12, "3D x stride = height * depth");
}
