#include "hydro/limiter/limiter.h"
#include "config.h"
#include "test_support.h"
#include <array>
#include <fstream>
#include <limits>

int main(int argc, char** argv)
{
    return run_test("limiter", [&] {
        CHECK(argc == 2);
        using hydro::Limiter;
        using hydro::limited_slope;
        CHECK(limited_slope(1, 3, Limiter::none) == 2);
        CHECK(limited_slope(1, -3, Limiter::none) == -1);
        CHECK(limited_slope(0, 3, Limiter::none) == 1.5);
        CHECK(limited_slope(-3, -1, Limiter::none) == -2);
        const double largest = std::numeric_limits<double>::max();
        CHECK(limited_slope(largest, largest, Limiter::none) == largest);
        expect_throw<std::invalid_argument>([] {
            limited_slope(std::numeric_limits<double>::infinity(), 1, Limiter::none);
        });
        const std::array<Limiter, 5> kinds{
            Limiter::minmod, Limiter::MC, Limiter::van_Leer,
            Limiter::MCbeta, Limiter::superbee
        };
        const std::array<double, 5> expected{1.0, 2.0, 1.5, 1.5, 2.0};
        for (std::size_t i = 0; i < kinds.size(); ++i)
        {
            const auto kind = kinds[i];
            CHECK(near(limited_slope(1.0, 3.0, kind), expected[i]));
            CHECK(near(limited_slope(3.0, 1.0, kind), expected[i]));
            CHECK(near(limited_slope(-1.0, -3.0, kind), -expected[i]));
            CHECK(limited_slope(1.0, -3.0, kind) == 0.0);
            CHECK(limited_slope(-1.0, 3.0, kind) == 0.0);
            CHECK(limited_slope(0.0, 3.0, kind) == 0.0);
            // Linear profiles must be reproduced, including extreme magnitudes.
            for (double slope : {0.0, 2.0, -2.0, 1e-250, 1e250,
                                 std::numeric_limits<double>::max()})
                CHECK(limited_slope(slope, slope, kind) == slope);
            const double extreme = limited_slope(1e250, 1.5e250, kind);
            CHECK(std::isfinite(extreme) && extreme >= 1e250 && extreme <= 1.5e250);
            expect_throw<std::invalid_argument>([&] {
                limited_slope(std::numeric_limits<double>::infinity(), 1.0, kind);
            });
            expect_throw<std::invalid_argument>([&] {
                limited_slope(1.0, std::numeric_limits<double>::quiet_NaN(), kind);
            });
        }
        // MC-beta endpoints recover minmod and MC for unequal slopes.
        CHECK(limited_slope(1, 3, Limiter::MCbeta, 1) == limited_slope(1, 3, Limiter::minmod));
        CHECK(limited_slope(1, 3, Limiter::MCbeta, 2) == limited_slope(1, 3, Limiter::MC));
        for (double beta : {0.9, 2.1, std::numeric_limits<double>::quiet_NaN()})
            expect_throw<std::invalid_argument>([&] { limited_slope(0, 0, Limiter::MCbeta, beta); });
        expect_throw<std::invalid_argument>([] { limited_slope(0, 0, static_cast<Limiter>(99)); });

        CHECK(hydro::parse_limiter("NONE") == Limiter::none);
        CHECK(hydro::parse_limiter("MINMOD") == Limiter::minmod);
        CHECK(hydro::parse_limiter("MC") == Limiter::MC);
        CHECK(hydro::parse_limiter("vanleer") == Limiter::van_Leer);
        CHECK(hydro::parse_limiter("van_Leer") == Limiter::van_Leer);
        CHECK(hydro::parse_limiter("mcbeta") == Limiter::MCbeta);
        CHECK(hydro::parse_limiter("mc_beta") == Limiter::MCbeta);
        CHECK(hydro::parse_limiter("superbee") == Limiter::superbee);
        expect_throw<std::invalid_argument>([] { hydro::parse_limiter("unknown"); });

        const std::string path = std::string(argv[1]) + ".txt";
        const auto read_config = [&](const std::string& extra) {
            {
                std::ofstream file(path);
                const char* axes[] = {"x", "y", "z"};
                for (int d = 0; d < hd_core::NDIM; ++d)
                    file << "n" << axes[d] << "a = 8\nn" << axes[d] << "g = 2\n"
                         << axes[d] << "_lower = 0\n" << axes[d] << "_upper = 1\n";
                file << extra;
                CHECK(file.good());
            }
            return config::ReadParameters(path);
        };
        CHECK(read_config("").limiter == "minmod");
        for (const std::string name : {"none", "minmod", "MC", "vanleer", "mcbeta", "superbee"})
            CHECK(read_config("limiter = " + name + "\n").limiter == name);
        const auto cfg = read_config("limiter = mcbeta\nlimiter_beta = 1.8\n");
        CHECK(cfg.limiter_beta == 1.8);
        CHECK(near(limited_slope(1, 3, hydro::parse_limiter(cfg.limiter), cfg.limiter_beta), 1.8));
        expect_throw<std::invalid_argument>([&] { read_config("limiter = typo\n"); });
        expect_throw<std::runtime_error>([&] { read_config("limiter_beta = 2.1\n"); });
    });
}
