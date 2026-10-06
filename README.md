# bromath

[![CI](https://github.com/wlejon/bromath/actions/workflows/ci.yml/badge.svg)](https://github.com/wlejon/bromath/actions/workflows/ci.yml)
[![CodeQL](https://github.com/wlejon/bromath/actions/workflows/codeql.yml/badge.svg)](https://github.com/wlejon/bromath/actions/workflows/codeql.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Shared math primitives for the bro stack. Header-only, C++20, no
third-party dependencies.

bromath is the bottom of the [bro ecosystem](https://github.com/wlejon/bro/blob/main/docs/ecosystem.md):
it depends on nothing, and these repos build against it directly: bro,
broaudio, bromesh, broflora, brogameagent, broimage, brolm, brodiffusion,
brosoundml, brovisionml, bromux and brothumb. Others (brokit, through
broimage) get it transitively. It has no JavaScript binding of its own; bro
exposes the parts apps need under `bro.math`.

Being header-only, it runs wherever a C++20 compiler does. `simd.h` uses
SSE2/AVX2 on x86-64 and NEON on AArch64, with a scalar path everywhere else.

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
| `easing.h` | Penner easing set: easeLinear and the In/Out/InOut variants of Quad, Cubic, Quart, Quint, Sine, Expo, Circ, Back, Elastic, Bounce |
| `simd.h` | Simd4f (SSE2/AVX2, NEON, scalar fallback; `BROMATH_SIMD_FORCE_SCALAR` forces scalar) and batch Vec3 kernels: batchAdd/Sub/Scale/MulAdd/Lerp/Dot/Dist2/Normalize/Min/Max, batchTransformPoints/Vectors, batchSphereOverlap, batchAABBContains |
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
- Headers `#include` only the C++ standard library (plus the platform
  SIMD intrinsics headers in `simd.h`). `color.h`'s CSS colour parser is
  the one user of `<string>` / `<string_view>`.

## Build

```bash
cmake -B build
cmake --build build --config Debug
ctest --test-dir build -C Debug --output-on-failure
```

CI builds and runs the tests on Linux (GCC + Clang), Windows (MSVC) and macOS/arm64. The matrix is the point of it: a header-only library is compiled fresh inside every consumer, under whatever toolchain that consumer uses, so a construct one compiler accepts and another rejects has to be caught here rather than in whichever sibling next builds on the other one.

Coverage of `include/bromath/` is reported in each run's job summary (`-DBROMATH_COVERAGE=ON` locally; GCC/Clang only). [CodeQL](.github/workflows/codeql.yml) analyses the headers weekly and on every push.

## Consuming bromath

Header-only INTERFACE library (`bromath`, alias `bromath::bromath`).
Consumers resolve it the way every repo in the ecosystem resolves a
sibling: an existing `bromath` target wins (a superbuild already added it),
then a checkout beside the top-level project at `../bromath`, then the
top-level project's `third_party/bromath` git submodule:

```cmake
set(BROMATH_DIR "${CMAKE_SOURCE_DIR}/../bromath" CACHE PATH "Standalone bromath checkout")
if(NOT TARGET bromath)
    if(EXISTS "${BROMATH_DIR}/CMakeLists.txt")
        add_subdirectory("${BROMATH_DIR}" "${CMAKE_BINARY_DIR}/bromath" EXCLUDE_FROM_ALL)
    else()
        add_subdirectory("${CMAKE_SOURCE_DIR}/third_party/bromath"
                         "${CMAKE_BINARY_DIR}/bromath" EXCLUDE_FROM_ALL)
    endif()
endif()

target_link_libraries(your_target PUBLIC bromath::bromath)
```

Tests are built only when bromath is the top-level project
(`BROMATH_TESTS` defaults to `PROJECT_IS_TOP_LEVEL`).

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

- **Tensor / NN math** — brotensor (unified Tensor type; CPU, CUDA, Metal and Vulkan ops)
- **DSP** (biquad, FFT, polyBLEP, resampler) — broaudio
- **Mesh operations** (CSG, remesh, simplify, raycast acceleration) — bromesh
- **Procedural noise** (Simplex, FBm) — FastNoise2, built by brokit and bro
- **Steering / AI** (seek/arrive/flee/pursue, intercept solver) — brogameagent
- **Spatial accel structures** (BVH) — bromesh

These may be extracted later if a second consumer materializes.

## Versioning

Pre-1.0. Siblings vendor this repo via `add_subdirectory` and compile the headers themselves, so a tag is a pin point for `FetchContent GIT_TAG` rather than a compatibility promise.

## License

[MIT](LICENSE)
