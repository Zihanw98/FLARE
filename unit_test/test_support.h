#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

// Throw instead of relying on NDEBUG-sensitive assertions in new tests.
#define CHECK(expression) do { if (!(expression)) throw std::runtime_error( \
    std::string(__FILE__) + ":" + std::to_string(__LINE__) + ": " #expression); } while (false)

inline bool near(double a, double b, double tolerance = 1e-11)
{
    return std::isfinite(a) && std::isfinite(b)
        && std::abs(a - b) <= tolerance * std::max({1.0, std::abs(a), std::abs(b)});
}

template<class Exception, class Function>
void expect_throw(Function function)
{
    bool caught = false;
    try { function(); }
    catch (const Exception&) { caught = true; }
    CHECK(caught);
}

template<class Function>
int run_test(const char* name, Function function)
{
    try
    {
        function();
        std::cout << name << " passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << name << " FAILED: " << error.what() << '\n';
        return 1;
    }
}
