#pragma once

// Small hash helpers — non-cryptographic, deterministic. FNV-1a for byte
// sequences, integer mixers for scalar inputs, and positionToCell for
// spatial hash structures (game grids, GridFootprint2D).

#include "bromath/vec.h"

#include <cmath>
#include <cstddef>
#include <cstdint>

namespace bromath {

// FNV-1a 32-bit.
inline uint32_t fnv1a32(const void* data, std::size_t len, uint32_t seed = 2166136261u) {
    const uint8_t* p = (const uint8_t*)data;
    uint32_t h = seed;
    for (std::size_t i = 0; i < len; ++i) {
        h ^= p[i];
        h *= 16777619u;
    }
    return h;
}

// 32-bit integer hash (Wang).
inline constexpr uint32_t hashU32(uint32_t x) {
    x = (x ^ 61u) ^ (x >> 16);
    x = x + (x << 3);
    x = x ^ (x >> 4);
    x = x * 0x27d4eb2du;
    x = x ^ (x >> 15);
    return x;
}

// 64-bit integer hash (splitmix64 finalizer, pure function form).
inline constexpr uint64_t hashU64(uint64_t x) {
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

// Combine two hashes (boost-style).
inline constexpr uint32_t hashCombine(uint32_t a, uint32_t b) {
    return a ^ (b + 0x9e3779b9u + (a << 6) + (a >> 2));
}

// Hash a 3D integer cell coordinate to a 32-bit value — uniform over the
// spatial-hash use case.
inline constexpr uint32_t cellHash(int32_t x, int32_t y, int32_t z) {
    uint32_t h = hashU32((uint32_t)x * 0x9e3779b1u);
    h = hashCombine(h, hashU32((uint32_t)y * 0x85ebca77u));
    h = hashCombine(h, hashU32((uint32_t)z * 0xc2b2ae3du));
    return h;
}

inline constexpr uint32_t cellHash(int32_t x, int32_t z) {
    uint32_t h = hashU32((uint32_t)x * 0x9e3779b1u);
    return hashCombine(h, hashU32((uint32_t)z * 0x85ebca77u));
}

// Integer cell of a coordinate, clamped to +-2^30 so a huge, infinite or NaN
// value (NaN goes to the low end) is not an undefined float-to-int cast.
inline int32_t cellCoordOf(float v, float cellSize) {
    constexpr float kLimit = 1073741824.0f;  // 2^30
    const float f = std::floor(v / cellSize);
    if (!(f > -kLimit)) return -(1 << 30);
    if (f > kLimit) return 1 << 30;
    return (int32_t)f;
}

// Hash a position into a cell index, useful for sparse spatial-hash tables.
// A bucketCount of 0 gives 0.
inline uint32_t positionToCell(Vec3 p, float cellSize, uint32_t bucketCount) {
    if (bucketCount == 0) return 0;
    return cellHash(cellCoordOf(p.x, cellSize), cellCoordOf(p.y, cellSize),
                    cellCoordOf(p.z, cellSize)) % bucketCount;
}

} // namespace bromath
