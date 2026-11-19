#include "test_support.h"
#include "utils/hdmath.h"
#include <limits>

namespace
{
double quadratic(double x) { return x * x - 2.0; }
double quadratic_derivative(double x) { return 2.0 * x; }
}

int main()
{
    return run_test("generic Newton-Raphson", [] {
        using hdmath::newton_raphson;
        CHECK(near(newton_raphson(&quadratic, &quadratic_derivative, 1.0),
                   std::sqrt(2.0)));

        // Captured coefficients: inversion of thermal + radiation energy density.
        // Normalize by the target energy so tolerance has the same meaning in
        // gas-dominated and radiation-dominated states and across energy units.
        constexpr double a = 7.5657e-15;
        for (const double C : {1e-4, 1e4})
        for (const double expected : {1e4, 1e6})
        for (const double scale : {1e-100, 1.0, 1e100})
        {
            const double q = C * expected + a * std::pow(expected, 4);
            auto energy = [=](double T) {
                return (scale * C * T + scale * a * std::pow(T, 4) - scale * q)
                       / (scale * q);
            };
            auto derivative = [=](double T) {
                return (scale * C + 4.0 * scale * a * std::pow(T, 3)) / (scale * q);
            };
            const double T = newton_raphson(energy, derivative, 2.0 * expected, 1e-12);
            CHECK(near(T, expected, 1e-12));
            CHECK(std::abs(energy(T)) <= 1e-12);
        }

        auto linear = [](double x) { return x - 3.0; };
        auto one = [](double) { return 1.0; };
        CHECK(newton_raphson(linear, one, -2.0, 1e-12, 1) == 3.0);
        // An initially converged value needs no derivative, even at a multiple root.
        CHECK(newton_raphson([](double x) { return x*x; },
            [](double) -> double { throw std::logic_error("must not evaluate"); }, 0.0) == 0.0);

        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double inf = std::numeric_limits<double>::infinity();
        for (const double bad : {nan, inf})
        {
            expect_throw<std::invalid_argument>([&] { newton_raphson(linear, one, bad); });
            expect_throw<std::runtime_error>([&] {
                newton_raphson([=](double) { return bad; }, one, 1.0);
            });
            expect_throw<std::runtime_error>([&] {
                newton_raphson(linear, [=](double) { return bad; }, 1.0);
            });
        }
        for (const double bad : {0.0, -1.0, nan, inf})
            expect_throw<std::invalid_argument>([&] { newton_raphson(linear, one, 1.0, bad); });
        for (const int bad : {0, -1})
            expect_throw<std::invalid_argument>([&] { newton_raphson(linear, one, 1.0, 1e-10, bad); });

        expect_throw<std::runtime_error>([&] {
            newton_raphson(quadratic, quadratic_derivative, 0.0); // Zero slope, nonzero residual.
        });
        expect_throw<std::runtime_error>([&] {
            newton_raphson(quadratic, quadratic_derivative, 1.0, 1e-12, 1);
        });
        expect_throw<std::runtime_error>([&] {
            newton_raphson([](double x) { return x*x*x - 2*x + 2; },
                           [](double x) { return 3*x*x - 2; }, 0.0, 1e-10, 10); // Cycle 0 -> 1 -> 0.
        });
        expect_throw<std::runtime_error>([&] {
            newton_raphson(one, [](double) { return 1e300; }, 1.0); // Tiny update is not convergence.
        });
        expect_throw<std::runtime_error>([&] {
            newton_raphson([](double) { return 1e300; },
                           [](double) { return 1e-300; }, 1.0); // Division overflow.
        });
        expect_throw<std::runtime_error>([&] {
            newton_raphson([](double) { return -std::numeric_limits<double>::max(); },
                           one, std::numeric_limits<double>::max()); // Addition overflow.
        });
        expect_throw<std::domain_error>([&] {
            newton_raphson([](double) -> double { throw std::domain_error("outside EoS domain"); },
                           one, 1.0);
        });
    });
}
