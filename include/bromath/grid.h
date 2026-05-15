#pragma once

// Spatial grid metadata. GridFootprint2D (lifted from broflora) describes a
// row-major 2D heightfield or occupancy grid sitting in world space. The
// helpers here let callers convert between world coordinates and cell
// indices without each subsystem rewriting the same math.

#include "bromath/vec.h"

#include <cmath>
#include <cstdint>

namespace bromath {

struct GridFootprint2D {
    Vec2 origin = {0, 0};
    float cellSize = 1.0f;
    int   width  = 0;  // cells along X
    int   depth  = 0;  // cells along Z (or Y for top-down 2D)
};

// Row-major index for a (col, row) cell. width is the X stride.
inline constexpr int gridIndex2D(const GridFootprint2D& g, int col, int row) {
    return row * g.width + col;
}

inline constexpr bool gridInBounds(const GridFootprint2D& g, int col, int row) {
    return col >= 0 && col < g.width && row >= 0 && row < g.depth;
}

// Convert world XY to integer cell coordinate. May return out-of-range
// indices — pair with gridInBounds when reading from a grid array.
inline void gridCellOf(const GridFootprint2D& g, Vec2 p, int& col, int& row) {
    float fx = (p.x - g.origin.x) / g.cellSize;
    float fy = (p.y - g.origin.y) / g.cellSize;
    col = (int)std::floor(fx);
    row = (int)std::floor(fy);
}

// Centre point of cell (col, row) in world coordinates.
inline Vec2 gridCellCenter(const GridFootprint2D& g, int col, int row) {
    return { g.origin.x + (col + 0.5f) * g.cellSize,
             g.origin.y + (row + 0.5f) * g.cellSize };
}

// Row-major index helper for a 3D grid stored as (x * depth * height + y * depth + z).
// Chosen to match broflora's ShadowGrid layout (qg[x*h+y*d+z] with x outer).
inline constexpr int gridIndex3D(int x, int y, int z, int depth, int height) {
    return (x * height + y) * depth + z;
}

} // namespace bromath
