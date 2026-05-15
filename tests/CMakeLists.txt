add_executable(bromath_test
    test_main.cpp
    test_scalar.cpp
    test_vec.cpp
    test_quat.cpp
    test_mat.cpp
    test_geom.cpp
)

target_link_libraries(bromath_test PRIVATE bromath)
target_compile_features(bromath_test PRIVATE cxx_std_20)

if(MSVC)
    target_compile_options(bromath_test PRIVATE /W4)
else()
    target_compile_options(bromath_test PRIVATE -Wall -Wextra)
endif()

add_test(NAME bromath_test COMMAND bromath_test)
