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
// High-performance compact layout:
// - All entries reside in a contiguous std::vector<Entry> for maximum cache locality.
// - Cells are managed via an open-addressing flat hash table with 64-bit spatial
//   bit-scrambled hashing (SplitMix64) and backward-shift deletion.
// - An entity ID reverse map provides strictly O(1) item removal without scanning cells.

#include "bromath/aabb.h"
#include "bromath/hash.h"
#include "bromath/sphere.h"
#include "bromath/vec.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
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
        clear();
    }

    void clear() {
        entries_.clear();
        std::fill(cellTable_.begin(), cellTable_.end(), CellSlot{});
        std::fill(idTable_.begin(), idTable_.end(), IdSlot{});
        occupiedCells_ = 0;
        occupiedIds_ = 0;
        maxRadius_ = 0.0f;
    }

    // Insert a point. Equivalent to inserting a zero-radius sphere.
    void insert(Vec3 p, int32_t id) {
        insertEntry(p, 0.0f, id);
    }

    // Insert a sphere. The grid records the maximum radius ever inserted; all
    // queries dilate by that amount so a large sphere whose center sits
    // outside the query volume is still found.
    void insert(Sphere s, int32_t id) {
        insertEntry(s.center, s.radius, id);
        if (s.radius > maxRadius_) maxRadius_ = s.radius;
    }

    // Remove every entry whose id matches in O(1) time via entity reverse map.
    // maxRadius_ is left as an upper bound (never tightened) — a stale
    // maxRadius only inflates query cost, never correctness.
    void remove(int32_t id) {
        uint32_t idx = findId(id);
        if (idx == INVALID_INDEX) return;

        // 1. Unlink idx from its cell list
        uint64_t key = entries_[idx].cellKey;
        size_t cSlot = findCellSlot(key);
        uint32_t prev = entries_[idx].prevInCell;
        uint32_t next = entries_[idx].nextInCell;
        if (prev != INVALID_INDEX) {
            entries_[prev].nextInCell = next;
        } else if (cSlot != SIZE_MAX) {
            cellTable_[cSlot].head = next;
        }
        if (next != INVALID_INDEX) {
            entries_[next].prevInCell = prev;
        }
        if (cSlot != SIZE_MAX && cellTable_[cSlot].head == INVALID_INDEX) {
            eraseCellSlot(cSlot);
        }

        // 2. Swap-and-pop with back entry if idx is not the last entry
        uint32_t backIdx = static_cast<uint32_t>(entries_.size() - 1);
        if (idx != backIdx) {
            entries_[idx] = entries_[backIdx];
            if (entries_[idx].prevInCell != INVALID_INDEX) {
                entries_[entries_[idx].prevInCell].nextInCell = idx;
            } else {
                size_t swappedCSlot = findCellSlot(entries_[idx].cellKey);
                if (swappedCSlot != SIZE_MAX) {
                    cellTable_[swappedCSlot].head = idx;
                }
            }
            if (entries_[idx].nextInCell != INVALID_INDEX) {
                entries_[entries_[idx].nextInCell].prevInCell = idx;
            }
            setId(entries_[idx].id, idx);
        }

        // 3. Remove target id from reverse map and pop back
        eraseId(id);
        entries_.pop_back();
    }

    // Append ids of entries whose extent lies within `radius` of `center`.
    // Existing contents of `out` are not cleared.
    void radiusQuery(Vec3 center, float radius, std::vector<int32_t>& out) const {
        if (radius <= 0.0f && maxRadius_ <= 0.0f) return;
        if (entries_.empty() || cellTable_.empty()) return;
        const float searchR = radius + maxRadius_;
        const int extent = static_cast<int>(std::ceil(searchR * invCell_));
        int cx, cy, cz;
        cellOf(center, cx, cy, cz);
        for (int dz = -extent; dz <= extent; ++dz) {
            for (int dy = -extent; dy <= extent; ++dy) {
                for (int dx = -extent; dx <= extent; ++dx) {
                    uint64_t key = makeKey(cx + dx, cy + dy, cz + dz);
                    size_t slot = findCellSlot(key);
                    if (slot == SIZE_MAX) continue;
                    uint32_t cur = cellTable_[slot].head;
                    while (cur != INVALID_INDEX) {
                        const Entry& e = entries_[cur];
                        const float reach = radius + e.r;
                        if (vdist2(e.p, center) <= reach * reach) {
                            out.push_back(e.id);
                        }
                        cur = e.nextInCell;
                    }
                }
            }
        }
    }

    // Append ids of entries that touch `box`. Point entries match when their
    // center lies in the box; sphere entries match when their sphere
    // intersects the box.
    void queryAABB(const AABB3& box, std::vector<int32_t>& out) const {
        if (aisEmpty(box) || entries_.empty() || cellTable_.empty()) return;
        const float pad = maxRadius_;
        int ix0, iy0, iz0, ix1, iy1, iz1;
        cellOf({box.min.x - pad, box.min.y - pad, box.min.z - pad}, ix0, iy0, iz0);
        cellOf({box.max.x + pad, box.max.y + pad, box.max.z + pad}, ix1, iy1, iz1);
        for (int iz = iz0; iz <= iz1; ++iz) {
            for (int iy = iy0; iy <= iy1; ++iy) {
                for (int ix = ix0; ix <= ix1; ++ix) {
                    uint64_t key = makeKey(ix, iy, iz);
                    size_t slot = findCellSlot(key);
                    if (slot == SIZE_MAX) continue;
                    uint32_t cur = cellTable_[slot].head;
                    while (cur != INVALID_INDEX) {
                        const Entry& e = entries_[cur];
                        if (e.r <= 0.0f) {
                            if (acontains(box, e.p)) out.push_back(e.id);
                        } else if (sphereTouchesAABB(e.p, e.r, box)) {
                            out.push_back(e.id);
                        }
                        cur = e.nextInCell;
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
        if (maxRadius <= 0.0f || entries_.empty() || cellTable_.empty()) return -1;
        const int extent = static_cast<int>(std::ceil(maxRadius * invCell_));
        int cx, cy, cz;
        cellOf(center, cx, cy, cz);
        int32_t bestId = -1;
        float bestD2 = maxRadius * maxRadius;
        for (int dz = -extent; dz <= extent; ++dz) {
            for (int dy = -extent; dy <= extent; ++dy) {
                for (int dx = -extent; dx <= extent; ++dx) {
                    uint64_t key = makeKey(cx + dx, cy + dy, cz + dz);
                    size_t slot = findCellSlot(key);
                    if (slot == SIZE_MAX) continue;
                    uint32_t cur = cellTable_[slot].head;
                    while (cur != INVALID_INDEX) {
                        const Entry& e = entries_[cur];
                        float d2 = vdist2(e.p, center);
                        if (d2 < bestD2) {
                            bestD2 = d2;
                            bestId = e.id;
                        }
                        cur = e.nextInCell;
                    }
                }
            }
        }
        return bestId;
    }

    size_t size() const { return entries_.size(); }
    float cellSize() const { return cellSize_; }
    float maxRadius() const { return maxRadius_; }

private:
    static constexpr uint32_t INVALID_INDEX = 0xFFFFFFFFu;
    static constexpr uint64_t EMPTY_CELL_KEY = 0xFFFFFFFFFFFFFFFFULL;

    struct Entry {
        Vec3 p;
        float r;
        int32_t id;
        uint32_t nextInCell;
        uint32_t prevInCell;
        uint64_t cellKey;
    };

    struct CellSlot {
        uint64_t key = EMPTY_CELL_KEY;
        uint32_t head = INVALID_INDEX;
    };

    struct IdSlot {
        int32_t id = 0;
        uint32_t entryIdx = INVALID_INDEX;
    };

    // Pack 3 x 21-bit signed ints into one 64-bit key (range ~+-1M cells).
    static constexpr uint64_t makeKey(int ix, int iy, int iz) {
        constexpr uint64_t mask = (1ULL << 21) - 1ULL;
        uint64_t a = static_cast<uint64_t>(static_cast<int64_t>(ix)) & mask;
        uint64_t b = static_cast<uint64_t>(static_cast<int64_t>(iy)) & mask;
        uint64_t c = static_cast<uint64_t>(static_cast<int64_t>(iz)) & mask;
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

    void insertEntry(Vec3 p, float r, int32_t id) {
        // If id already exists, replace it cleanly
        if (findId(id) != INVALID_INDEX) {
            remove(id);
        }

        int ix, iy, iz;
        cellOf(p, ix, iy, iz);
        uint64_t key = makeKey(ix, iy, iz);

        size_t cSlot = findOrCreateCellSlot(key);
        uint32_t newIdx = static_cast<uint32_t>(entries_.size());
        uint32_t oldHead = cellTable_[cSlot].head;

        entries_.push_back({p, r, id, oldHead, INVALID_INDEX, key});
        if (oldHead != INVALID_INDEX) {
            entries_[oldHead].prevInCell = newIdx;
        }
        cellTable_[cSlot].head = newIdx;

        setId(id, newIdx);
    }

    // Cell flat hash table operations (linear probing + SplitMix64)
    size_t findCellSlot(uint64_t key) const {
        if (cellTable_.empty()) return SIZE_MAX;
        size_t mask = cellTable_.size() - 1;
        size_t i = hashU64(key) & mask;
        while (cellTable_[i].key != EMPTY_CELL_KEY) {
            if (cellTable_[i].key == key) return i;
            i = (i + 1) & mask;
        }
        return SIZE_MAX;
    }

    size_t findOrCreateCellSlot(uint64_t key) {
        if (cellTable_.empty() || (occupiedCells_ + 1) * 2 >= cellTable_.size()) {
            rehashCells(cellTable_.empty() ? 64 : cellTable_.size() * 2);
        }
        size_t mask = cellTable_.size() - 1;
        size_t i = hashU64(key) & mask;
        while (cellTable_[i].key != EMPTY_CELL_KEY) {
            if (cellTable_[i].key == key) return i;
            i = (i + 1) & mask;
        }
        cellTable_[i].key = key;
        cellTable_[i].head = INVALID_INDEX;
        ++occupiedCells_;
        return i;
    }

    void eraseCellSlot(size_t i) {
        if (cellTable_.empty()) return;
        size_t mask = cellTable_.size() - 1;
        size_t j = i;
        while (true) {
            j = (j + 1) & mask;
            if (cellTable_[j].key == EMPTY_CELL_KEY) break;
            size_t k = hashU64(cellTable_[j].key) & mask;
            if ((j > i) ? (k <= i || k > j) : (k <= i && k > j)) {
                cellTable_[i] = cellTable_[j];
                i = j;
            }
        }
        cellTable_[i].key = EMPTY_CELL_KEY;
        cellTable_[i].head = INVALID_INDEX;
        --occupiedCells_;
    }

    void rehashCells(size_t newCap) {
        std::vector<CellSlot> oldTable = std::move(cellTable_);
        cellTable_.assign(newCap, CellSlot{});
        occupiedCells_ = 0;
        size_t mask = newCap - 1;
        for (const auto& slot : oldTable) {
            if (slot.key != EMPTY_CELL_KEY && slot.head != INVALID_INDEX) {
                size_t i = hashU64(slot.key) & mask;
                while (cellTable_[i].key != EMPTY_CELL_KEY) {
                    i = (i + 1) & mask;
                }
                cellTable_[i] = slot;
                ++occupiedCells_;
            }
        }
    }

    // Entity ID flat reverse map operations (linear probing + Wang hash)
    uint32_t findId(int32_t id) const {
        if (idTable_.empty()) return INVALID_INDEX;
        size_t mask = idTable_.size() - 1;
        size_t i = hashU32(static_cast<uint32_t>(id)) & mask;
        while (idTable_[i].entryIdx != INVALID_INDEX) {
            if (idTable_[i].id == id) return idTable_[i].entryIdx;
            i = (i + 1) & mask;
        }
        return INVALID_INDEX;
    }

    void setId(int32_t id, uint32_t entryIdx) {
        if (idTable_.empty() || (occupiedIds_ + 1) * 2 >= idTable_.size()) {
            rehashIds(idTable_.empty() ? 64 : idTable_.size() * 2);
        }
        size_t mask = idTable_.size() - 1;
        size_t i = hashU32(static_cast<uint32_t>(id)) & mask;
        while (idTable_[i].entryIdx != INVALID_INDEX) {
            if (idTable_[i].id == id) {
                idTable_[i].entryIdx = entryIdx;
                return;
            }
            i = (i + 1) & mask;
        }
        idTable_[i].id = id;
        idTable_[i].entryIdx = entryIdx;
        ++occupiedIds_;
    }

    void eraseId(int32_t id) {
        if (idTable_.empty()) return;
        size_t mask = idTable_.size() - 1;
        size_t i = hashU32(static_cast<uint32_t>(id)) & mask;
        while (idTable_[i].entryIdx != INVALID_INDEX) {
            if (idTable_[i].id == id) {
                size_t j = i;
                while (true) {
                    j = (j + 1) & mask;
                    if (idTable_[j].entryIdx == INVALID_INDEX) break;
                    size_t k = hashU32(static_cast<uint32_t>(idTable_[j].id)) & mask;
                    if ((j > i) ? (k <= i || k > j) : (k <= i && k > j)) {
                        idTable_[i] = idTable_[j];
                        i = j;
                    }
                }
                idTable_[i].entryIdx = INVALID_INDEX;
                --occupiedIds_;
                return;
            }
            i = (i + 1) & mask;
        }
    }

    void rehashIds(size_t newCap) {
        std::vector<IdSlot> oldTable = std::move(idTable_);
        idTable_.assign(newCap, IdSlot{});
        occupiedIds_ = 0;
        size_t mask = newCap - 1;
        for (const auto& slot : oldTable) {
            if (slot.entryIdx != INVALID_INDEX) {
                size_t i = hashU32(static_cast<uint32_t>(slot.id)) & mask;
                while (idTable_[i].entryIdx != INVALID_INDEX) {
                    i = (i + 1) & mask;
                }
                idTable_[i] = slot;
                ++occupiedIds_;
            }
        }
    }

    float cellSize_ = 1.0f;
    float invCell_ = 1.0f;
    float maxRadius_ = 0.0f;

    std::vector<Entry> entries_;
    std::vector<CellSlot> cellTable_;
    std::vector<IdSlot> idTable_;
    size_t occupiedCells_ = 0;
    size_t occupiedIds_ = 0;
};

} // namespace bromath
