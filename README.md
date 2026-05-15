# bromath

Shared math primitives for the bro stack (bro, bromesh, brogameagent,
broaudio, broflora). Header-only, C++20, no third-party dependencies.

## Scope

Geometric, scalar, and small-data math used across multiple sibling
libraries. **Not** a home for domain runtimes — NN/tensor math lives in
brogameagent, DSP in broaudio, mesh operations in bromesh.

| Header | Contents |
|--------|----------|
| `scalar.h` | constants (PI, TWO_PI, HALF_PI, DEG2RAD, ...), clamp, saturate, lerp, mix, smoothstep, smootherstep, step, sign, sqr, deg2rad, rad2deg |
| `angle.h` | wrapAngle, wrapAngle2Pi, angleDelta, angleLerp |
| `vec.h` | Vec2, Vec3 + free-function ops (vdot, vcross, vlen, vnorm, vlerp, vreflect, vproject, vperpendicular, ...) |
| `quat.h` | Quat (xyzw) + qmul, qrotate, qaxisAngle, qfromTo, qfromEuler, qtoEuler, qslerp, qnlerp |
| `mat.h` | Mat4 (column-major) + mmul, minverse, mtranspose, mfromTRS, mdecompose, mlookAt, mperspective, mortho |
| `transform.h` | Transform { pos, rot, scale } with ttoMat4 / tfromMat4 / tmul composition |
| `aabb.h` | AABB2, AABB3, contains, intersects, expand, merge, fromPoints, transform (Arvo) |
| `plane.h` | Plane (implicit form), pfromPointNormal, pfromPoints, psignedDistance, pproject |
| `sphere.h` | Sphere + sintersects + sintersectVolume (lens-volume closed form) |
| `ray.h` | Ray + ray-vs-AABB (slab), ray-vs-sphere, ray-vs-plane, ray-vs-triangle (Möller-Trumbore) |
| `frustum.h` | Frustum (six planes from VP matrix via Gribb-Hartmann), point/AABB/sphere culling |
| `color.h` | Color (linear RGBA float), Color8 (sRGB byte), HSV, hex parser, sRGB↔linear |
| `curves.h` | CubicEase (CSS-style), cbezier, chermite, ccatmullRom (centripetal) |
| `rng.h` | SplitMix64 + randFloat01, randSigned, randRange, randInt, randNormal, randInUnitDisc/Sphere, randOnUnitSphere |
| `hash.h` | FNV-1a, hashU32, hashU64, hashCombine, cellHash, positionToCell |
| `smoother.h` | One-pole parameter smoother |
| `grid.h` | GridFootprint2D + 2D/3D index helpers |

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
- All headers `#include` only `<cmath>`, `<cstdint>`, `<cstddef>`,
  `<vector>`, `<limits>` from the stdlib.

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

- **Tensor / NN math** — brogameagent (`nn/ops.h`, CUDA kernels)
- **DSP** (biquad, FFT, polyBLEP, resampler) — broaudio
- **Mesh operations** (CSG, remesh, simplify, raycast acceleration) — bromesh
- **Procedural noise** (Simplex, FBm) — FastNoise2, vendored in bromesh
- **Steering / AI** (seek/arrive/flee/pursue, intercept solver) — brogameagent
- **Spatial accel structures** (BVH, SpatialHash3D) — bromesh

These may be extracted later if a second consumer materializes.
