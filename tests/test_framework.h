#pragma once

// Shared test harness for bromath_test. Custom TEST(name) macro registers
// tests at static-init time; main() in test_main.cpp dispatches them.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

extern int tests_run;
extern int tests_passed;

using TestFn = void(*)();
struct TestEntry { const char* name; TestFn fn; };
std::vector<TestEntry>& testRegistry();

struct TestRegistrar {
    TestRegistrar(const char* name, TestFn fn) {
        testRegistry().push_back({name, fn});
    }
};

#define TEST(name) \
    static void test_##name(); \
    static TestRegistrar reg_##name(#name, &test_##name); \
    static void test_##name()

#define ASSERT(cond, msg) do { \
    tests_run++; \
    if (!(cond)) { \
        std::fprintf(stderr, "FAIL: %s (line %d)\n", msg, __LINE__); \
    } else { \
        tests_passed++; \
    } \
} while(0)

inline bool nearly(float a, float b, float eps = 1e-5f) {
    return std::fabs(a - b) <= eps;
}
