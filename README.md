# bromath

Shared math primitives for the bro stack — used transitively by most
siblings (bro, bromesh, brogameagent, broaudio, broflora, brotensor,
brolm, brodiffusion, broimage, brosoundml, brovisionml). Header-only,
C++20, no third-party dependencies.

## Scope

Geometric, scalar, and small-data math used across multiple sibling
libraries. **Not** a home for domain runtimes — NN/tensor math lives in
brotensor, DSP in broaudio, mesh operations in bromesh.

| Header | Contents |
|--------|----------|
| `scalar.h` | constants (PI, TWO_PI, HALF_PI, INV_PI, INV_TWO_PI, DEG2RAD, RAD2DEG), min, max, clamp, saturate, lerp, mix, invLerp, remap, smoothstep, smoothstep01, smootherstep, step, sign, abs, sqr, nearlyEqual, deg2rad, rad2deg |
| `angle.h` | wrapAngle, wrapAngle2Pi, angleDelta, angleLerp |
| `vec.h` | Vec2, Vec3 + free-function ops (vdot, vcross, vlen, vlen2, vnorm, vnormOr, vdist, vdist2, vlerp, vreflect, vproject, vperpendicular, vmin, vmax) |
| `quat.h` | Quat (xyzw) + qidentity, qmul, qconjugate, qinverse, qdot, qlen, qlen2, qnorm, qrotate, qaxisAngle, qfromTo, qfromEuler, qtoEuler, qslerp, qnlerp |
| `mat.h` | Mat4 (column-major) + midentity, mmul, minverse, mtranspose, mtranslate, mscale, mfromQuat, mfromTRS, mdecompose, mtransformPoint, mtransformDir, mlookAt, mperspective, mortho |
| `transform.h` | Transform { pos, rot, scale } with tidentity, ttoMat4, tfromMat4, tmul, ttransformPoint, ttransformDir |
| `aabb.h` | AABB2, AABB3 + acontains, aintersects, aexpand, amerge, afromPoints, atransform (Arvo), acenter, aextent, ahalfExtent, aisEmpty |
| `plane.h` | Plane (implicit form), pfromPointNormal, pfromPoints, psignedDistance, pproject |
| `sphere.h` | Sphere + scontains, sintersects, sintersectVolume (lens-volume closed form) |
| `segment.h` | Capsule (segment + radius) + closestSegmentSegment2, segmentSegmentDistance, capsulePenetration, capsulesIntersect (Ericson §5.1.9) |
| `ray.h` | Ray, RayHit + rat, rIntersectAABB (slab), rIntersectSphere, rIntersectPlane, rIntersectTriangle (Möller-Trumbore) |
| `frustum.h` | Frustum (six planes from VP matrix via Gribb-Hartmann) + ffromViewProj, fcontains, fintersects (point/AABB/sphere culling) |
| `color.h` | Color (linear RGBA float), Color8 (sRGB byte), cfromHSV, cfromHex, cfromColor8, ctoColor8, clerp, csrgbToLinear, clinearToSrgb |
| `curves.h` | ccubicEase (CSS-style), cbezier, cbezierTangent, chermite, ccatmullRom (centripetal) |
| `rng.h` | splitmix64 + randFloat01, randSigned, randRange, randInt, randNormal, randGaussian2D, randInUnitDisc, randInUnitSphere, randOnUnitSphere |
| `hash.h` | fnv1a32, hashU32, hashU64, hashCombine, cellHash, positionToCell |
| `smoother.h` | One-pole parameter smoother (smootherReset, smootherTarget, smootherSetTime, smootherTick, smootherTickN) |
| `grid.h` | GridFootprint2D + 2D/3D index helpers (gridIndex2D, gridIndex3D, gridCellOf, gridCellCenter, gridInBounds) |
| `spatial_hash.h` | SpatialHash3D (point and sphere indexing, radius and AABB queries) |
| `bromath.h` | umbrella header pulling in all of the above |

## Conventions

- **Free functions** for vector/matrix/quaternion ops: `vdot(a, b)`,
  `qrotate(q, v)`, `mmul(a, b)`. POD aggregates stay trivially copyable
  and bindings-friendly.
- **Matrices are column-major** (OpenGL / glTF / Jolt). `Mat4::data` is
  suitable for `glUniformMatrix4fv` with `transpose=GL_FALSE`.
- **Quaternions are xyzw** with identity `(0,0,0,1)`.
- **Vec2 is XY**. Other conventions (XZ for top-down nav) stay local to
  the consuming library.
- **Angles are radians** unless explicitly named otherwise.
- Headers `#include` only `<cmath>`, `<cstdint>`, `<limits>`,
  `<vector>`, and (spatial_hash only) `<unordered_map>` from the stdlib.

## Build

```bash
cmake -B build
cmake --build build --config Debug
./build/tests/Debug/bromath_test.exe
```

## Consuming bromath

Header-only INTERFACE library. From a sibling CMakeLists:

```cmake
# Prefer standalone repo, fall back to submodule
set(BROMATH_DIR "${CMAKE_SOURCE_DIR}/../bromath" CACHE PATH "")
if(EXISTS "${BROMATH_DIR}/CMakeLists.txt")
    add_subdirectory("${BROMATH_DIR}" "${CMAKE_BINARY_DIR}/bromath" EXCLUDE_FROM_ALL)
else()
    add_subdirectory(third_party/bromath EXCLUDE_FROM_ALL)
endif()

target_link_libraries(your_target PUBLIC bromath::bromath)
```

Then in code:

```cpp
#include <bromath/vec.h>
#include <bromath/quat.h>

using bromath::Vec3;
using bromath::Quat;
using bromath::vdot;
```

## Out of scope

The following intentionally live elsewhere:

- **Tensor / NN math** — brotensor (unified Tensor type, CPU/CUDA/Metal ops)
- **DSP** (biquad, FFT, polyBLEP, resampler) — broaudio
- **Mesh operations** (CSG, remesh, simplify, raycast acceleration) — bromesh
- **Procedural noise** (Simplex, FBm) — FastNoise2, vendored in bromesh
- **Steering / AI** (seek/arrive/flee/pursue, intercept solver) — brogameagent
- **Spatial accel structures** (BVH) — bromesh

These may be extracted later if a second consumer materializes.
