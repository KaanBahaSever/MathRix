#pragma once

// Minimal dependency-free test framework.

#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace mxtest {

struct Case {
    const char* name;
    std::function<void()> fn;
};

inline std::vector<Case>& registry() {
    static std::vector<Case> r;
    return r;
}

struct Registrar {
    Registrar(const char* name, std::function<void()> fn) { registry().push_back({name, std::move(fn)}); }
};

struct Failure {
    std::string message;
};

inline int& checks() {
    static int n = 0;
    return n;
}

inline void fail(const char* file, int line, const std::string& msg) {
    std::ostringstream ss;
    ss << file << ":" << line << ": " << msg;
    throw Failure{ss.str()};
}

inline int runAll(const std::string& filter = "") {
    int failed = 0, ran = 0;
    for (const auto& c : registry()) {
        if (!filter.empty() && std::string(c.name).find(filter) == std::string::npos) continue;
        ++ran;
        try {
            c.fn();
            std::cout << "[ OK ] " << c.name << "\n";
        } catch (const Failure& f) {
            ++failed;
            std::cout << "[FAIL] " << c.name << "\n       " << f.message << "\n";
        } catch (const std::exception& e) {
            ++failed;
            std::cout << "[FAIL] " << c.name << "\n       exception: " << e.what() << "\n";
        }
    }
    std::cout << "\n" << ran - failed << "/" << ran << " tests passed (" << checks() << " checks)\n";
    return failed == 0 ? 0 : 1;
}

}  // namespace mxtest

#define MX_CONCAT2(a, b) a##b
#define MX_CONCAT(a, b) MX_CONCAT2(a, b)
#define TEST_CASE(name)                                                                         \
    static void MX_CONCAT(mx_test_fn_, __LINE__)();                                             \
    static ::mxtest::Registrar MX_CONCAT(mx_test_reg_, __LINE__)(name, &MX_CONCAT(mx_test_fn_, __LINE__)); \
    static void MX_CONCAT(mx_test_fn_, __LINE__)()

#define CHECK(cond)                                                                       \
    do {                                                                                  \
        ++::mxtest::checks();                                                             \
        if (!(cond)) ::mxtest::fail(__FILE__, __LINE__, std::string("CHECK failed: ") + #cond); \
    } while (0)

#define CHECK_NEAR(a, b, tol)                                                                          \
    do {                                                                                               \
        ++::mxtest::checks();                                                                          \
        const double mx_a = (a), mx_b = (b);                                                           \
        if (!(std::abs(mx_a - mx_b) <= (tol))) {                                                       \
            std::ostringstream mx_ss;                                                                  \
            mx_ss << "CHECK_NEAR failed: " #a " = " << mx_a << ", " #b " = " << mx_b << ", tol " << (tol); \
            ::mxtest::fail(__FILE__, __LINE__, mx_ss.str());                                           \
        }                                                                                              \
    } while (0)

#define CHECK_REL(a, b, rel) CHECK_NEAR(a, b, std::abs(b) * (rel))

#define CHECK_THROWS(expr)                                                     \
    do {                                                                       \
        ++::mxtest::checks();                                                  \
        bool mx_thrown = false;                                                \
        try {                                                                  \
            (void)(expr);                                                      \
        } catch (...) {                                                        \
            mx_thrown = true;                                                  \
        }                                                                      \
        if (!mx_thrown) ::mxtest::fail(__FILE__, __LINE__, "expected exception: " #expr); \
    } while (0)
