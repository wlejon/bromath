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

# Header hygiene: every public header compiles on its own, as the first and
# only include of a translation unit, so none leans on what another brought in.
file(GLOB BROMATH_PUBLIC_HEADERS CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/include/bromath/*.h")
set(BROMATH_HEADER_TUS)
foreach(header IN LISTS BROMATH_PUBLIC_HEADERS)
    get_filename_component(name "${header}" NAME_WE)
    set(tu "${CMAKE_CURRENT_BINARY_DIR}/header_check/${name}.cpp")
    file(CONFIGURE OUTPUT "${tu}" CONTENT "#include \"bromath/${name}.h\"\n")
    list(APPEND BROMATH_HEADER_TUS "${tu}")
endforeach()
add_library(bromath_header_check OBJECT ${BROMATH_HEADER_TUS})
target_link_libraries(bromath_header_check PRIVATE bromath)
target_compile_features(bromath_header_check PRIVATE cxx_std_20)
if(MSVC)
    target_compile_options(bromath_header_check PRIVATE /W4)
else()
    target_compile_options(bromath_header_check PRIVATE -Wall -Wextra)
endif()
add_dependencies(bromath_test bromath_header_check)

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
