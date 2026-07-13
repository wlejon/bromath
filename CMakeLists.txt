cmake_minimum_required(VERSION 3.24)

# Pre-1.0: the version is a label for a tag, so consumers have something to pin
# a FetchContent GIT_TAG to. Siblings vendor this repo via add_subdirectory and
# compile the headers themselves, so there is no ABI to promise.
project(bromath VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Coverage. There is no library TU to instrument — the headers are the library —
# so this instruments the tests, and the headers get counted where the test TUs
# include them. Scoped to a top-level build: a consumer that add_subdirectory()s
# bromath has no business having its own build instrumented from in here.
option(BROMATH_COVERAGE "Instrument the tests for gcov coverage (GCC/Clang)" OFF)
if(BROMATH_COVERAGE AND PROJECT_IS_TOP_LEVEL)
    if(MSVC)
        message(FATAL_ERROR
            "BROMATH_COVERAGE requires GCC or Clang (gcov instrumentation); MSVC is unsupported.")
    endif()
    add_compile_options(--coverage -O0 -g)
    add_link_options(--coverage)
    message(STATUS "bromath: coverage instrumentation ON (-O0)")
endif()

if(MSVC)
    add_compile_options(/MP /permissive- /utf-8)
    if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()
endif()

# Header-only library. Consumers link bromath::bromath and #include <bromath/...>.
add_library(bromath INTERFACE)
add_library(bromath::bromath ALIAS bromath)

target_include_directories(bromath INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(bromath INTERFACE cxx_std_20)

option(BROMATH_TESTS "Build bromath tests" ${PROJECT_IS_TOP_LEVEL})
if(BROMATH_TESTS AND CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    enable_testing()
    add_subdirectory(tests)
endif()
