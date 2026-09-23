add_executable(bromath_test
    test_main.cpp
    test_scalar.cpp
    test_vec.cpp
    test_quat.cpp
    test_mat.cpp
    test_geom.cpp
    test_frustum.cpp
    test_color_curves.cpp
    test_easing.cpp
    test_rng_hash.cpp
    test_smoother_grid.cpp
    test_spatial_hash.cpp
    test_simd.cpp
    test_umbrella.cpp
)

target_link_libraries(bromath_test PRIVATE bromath)
target_compile_features(bromath_test PRIVATE cxx_std_20)

if(MSVC)
    target_compile_options(bromath_test PRIVATE /W4)
else()
    target_compile_options(bromath_test PRIVATE -Wall -Wextra)
endif()

add_test(NAME bromath_test COMMAND bromath_test)

# The SIMD tests again on the scalar fallback, so its results are checked
# against the same expectations as the intrinsic path.
add_executable(bromath_test_simd_scalar test_main.cpp test_simd.cpp)
target_link_libraries(bromath_test_simd_scalar PRIVATE bromath)
target_compile_features(bromath_test_simd_scalar PRIVATE cxx_std_20)
target_compile_definitions(bromath_test_simd_scalar PRIVATE BROMATH_SIMD_FORCE_SCALAR)
if(MSVC)
    target_compile_options(bromath_test_simd_scalar PRIVATE /W4)
else()
    target_compile_options(bromath_test_simd_scalar PRIVATE -Wall -Wextra)
endif()
add_test(NAME bromath_test_simd_scalar COMMAND bromath_test_simd_scalar)
