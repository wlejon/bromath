#pragma once

// Uniform-grid 3D spatial hash. Indexes points and spheres with int32_t ids;
// supports radius queries, AABB queries, and nearest-point lookup. Single-
// threaded by design — the project forbids mutexes; share by copying or by
// running independent grids per worker.
//
// Point vs sphere inserts: a point entry has radius 0 and matches a radius
// query when dist(p, center) <= radius. A sphere entry matches when its
// surface comes within `radius` of the query center, i.e.
// dist(p, center) <= radius + r. The grid tracks the maximum radius ever
// inserted and dilates the cell footprint of every query by that amount, so
// a large sphere whose center sits several cells outside the query volume
// is still considered.
//
// Memory: O(N) entries, plus one unordered_map<int64_t, vector<Entry>> bucket
// per occupied cell. Choose cellSize close to the typical query radius for
// best performance; very small cells inflate map churn, very large cells
// degenerate to linear scans.

#include "bromath/aabb.h"
#include "bromath/sphere.h"
#include "bromath/vec.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

namespace bromath {

class SpatialHash3D {
public:
    explicit SpatialHash3D(float cellSize = 1.0f) { reset(cellSize); }

    // Reset cell size and clear all entries. cellSize must be > 0; non-positive
    // values are clamped to 1.
    void reset(float cellSize) {
        cellSize_ = (cellSize > 0.0f) ? cellSize : 1.0f;
        invCell_ = 1.0f / cellSize_;
        cells_.clear();
        count_ = 0;
        maxRadius_ = 0.0f;
    }

    void clear() {
        cells_.clear();
        count_ = 0;
        maxRadius_ = 0.0f;
    }

    // Insert a point. Equivalent to inserting a zero-radius sphere.
    void insert(Vec3 p, int32_t id) {
        int ix, iy, iz;
        cellOf(p, ix, iy, iz);
        cells_[makeKey(ix, iy, iz)].push_back({p, 0.0f, id});
        ++count_;
    }

    // Insert a sphere. The grid records the maximum radius ever inserted; all
    // queries dilate by that amount so a large sphere whose center sits
    // outside the query volume is still found.
    void insert(Sphere s, int32_t id) {
        int ix, iy, iz;
        cellOf(s.center, ix, iy, iz);
        cells_[makeKey(ix, iy, iz)].push_back({s.center, s.radius, id});
        ++count_;
        if (s.radius > maxRadius_) maxRadius_ = s.radius;
    }

    // Remove every entry whose id matches. O(N) over touched cells.
    // maxRadius_ is left as an upper bound (never tightened) — a stale
    // maxRadius only inflates query cost, never correctness.
    void remove(int32_t id) {
        for (auto it = cells_.begin(); it != cells_.end(); ) {
            auto& vec = it->second;
            for (size_t i = vec.size(); i > 0; --i) {
                if (vec[i - 1].id == id) {
                    vec[i - 1] = vec.back();
                    vec.pop_back();
                    --count_;
                }
            }
            if (vec.empty()) it = cells_.erase(it);
            else ++it;
        }
    }

    // Append ids of entries whose extent lies within `radius` of `center`.
    // Existing contents of `out` are not cleared.
    void radiusQuery(Vec3 center, float radius, std::vector<int32_t>& out) const {
        if (radius <= 0.0f && maxRadius_ <= 0.0f) return;
        const float searchR = radius + maxRadius_;
        const int extent = static_cast<int>(std::ceil(searchR * invCell_));
        int cx, cy, cz;
        cellOf(center, cx, cy, cz);
        for (int dz = -extent; dz <= extent; ++dz) {
            for (int dy = -extent; dy <= extent; ++dy) {
                for (int dx = -extent; dx <= extent; ++dx) {
                    auto it = cells_.find(makeKey(cx + dx, cy + dy, cz + dz));
                    if (it == cells_.end()) continue;
                    for (const Entry& e : it->second) {
                        const float reach = radius + e.r;
                        if (vdist2(e.p, center) <= reach * reach) out.push_back(e.id);
                    }
                }
            }
        }
    }

    // Append ids of entries that touch `box`. Point entries match when their
    // center lies in the box; sphere entries match when their sphere
    // intersects the box.
    void queryAABB(const AABB3& box, std::vector<int32_t>& out) const {
        if (aisEmpty(box)) return;
        const float pad = maxRadius_;
        int ix0, iy0, iz0, ix1, iy1, iz1;
        cellOf({box.min.x - pad, box.min.y - pad, box.min.z - pad}, ix0, iy0, iz0);
        cellOf({box.max.x + pad, box.max.y + pad, box.max.z + pad}, ix1, iy1, iz1);
        for (int iz = iz0; iz <= iz1; ++iz) {
            for (int iy = iy0; iy <= iy1; ++iy) {
                for (int ix = ix0; ix <= ix1; ++ix) {
                    auto it = cells_.find(makeKey(ix, iy, iz));
                    if (it == cells_.end()) continue;
                    for (const Entry& e : it->second) {
                        if (e.r <= 0.0f) {
                            if (acontains(box, e.p)) out.push_back(e.id);
                        } else if (sphereTouchesAABB(e.p, e.r, box)) {
                            out.push_back(e.id);
                        }
                    }
                }
            }
        }
    }

    // Nearest entry whose center lies within `maxRadius` of `center`. Uses
    // center-to-center distance only — sphere radii are ignored. For sphere
    // semantics, use radiusQuery and pick the smallest reach yourself.
    // Returns -1 if no entry is within range.
    int32_t nearest(Vec3 center, float maxRadius) const {
        if (maxRadius <= 0.0f) return -1;
        const int extent = static_cast<int>(std::ceil(maxRadius * invCell_));
        int cx, cy, cz;
        cellOf(center, cx, cy, cz);
        int32_t bestId = -1;
        float bestD2 = maxRadius * maxRadius;
        for (int dz = -extent; dz <= extent; ++dz) {
            for (int dy = -extent; dy <= extent; ++dy) {
                for (int dx = -extent; dx <= extent; ++dx) {
                    auto it = cells_.find(makeKey(cx + dx, cy + dy, cz + dz));
                    if (it == cells_.end()) continue;
                    for (const Entry& e : it->second) {
                        float d2 = vdist2(e.p, center);
                        if (d2 < bestD2) {
                            bestD2 = d2;
                            bestId = e.id;
                        }
                    }
                }
            }
        }
        return bestId;
    }

    size_t size() const { return count_; }
    float cellSize() const { return cellSize_; }
    float maxRadius() const { return maxRadius_; }

private:
    struct Entry {
        Vec3 p;
        float r;
        int32_t id;
    };
    using CellKey = int64_t;

    // Pack 3 x 21-bit signed ints into one 64-bit key (range ~+-1M cells).
    static CellKey makeKey(int ix, int iy, int iz) {
        constexpr int64_t mask = (1LL << 21) - 1;
        int64_t a = static_cast<int64_t>(ix) & mask;
        int64_t b = static_cast<int64_t>(iy) & mask;
        int64_t c = static_cast<int64_t>(iz) & mask;
        return a | (b << 21) | (c << 42);
    }

    void cellOf(Vec3 p, int& ix, int& iy, int& iz) const {
        ix = static_cast<int>(std::floor(p.x * invCell_));
        iy = static_cast<int>(std::floor(p.y * invCell_));
        iz = static_cast<int>(std::floor(p.z * invCell_));
    }

    static bool sphereTouchesAABB(Vec3 c, float r, const AABB3& b) {
        float dx = c.x < b.min.x ? b.min.x - c.x : (c.x > b.max.x ? c.x - b.max.x : 0.0f);
        float dy = c.y < b.min.y ? b.min.y - c.y : (c.y > b.max.y ? c.y - b.max.y : 0.0f);
        float dz = c.z < b.min.z ? b.min.z - c.z : (c.z > b.max.z ? c.z - b.max.z : 0.0f);
        return dx * dx + dy * dy + dz * dz <= r * r;
    }

    float cellSize_ = 1.0f;
    float invCell_ = 1.0f;
    float maxRadius_ = 0.0f;
    size_t count_ = 0;
    std::unordered_map<CellKey, std::vector<Entry>> cells_;
};

} // namespace bromath
