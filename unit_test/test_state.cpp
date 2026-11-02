#include "test_support.h"
#include "core/cell.h"
#include "hydro/flux.h"
#include <limits>

int main()
{
    return run_test("state conversion, arithmetic and physical flux", [] {
        constexpr double gamma = 1.4;
        hd_core::Primitives w{};
        w[hd_core::RHO] = 2.0;
        w[hd_core::PRE] = 3.0;
        double speed_squared = 0.0;
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            const double v = (d % 2 == 0 ? 1 : -1) * (d + 1.0);
            w[hd_core::VEL(d)] = v;
            speed_squared += v * v;
        }
        const auto u = w.to_Conserved(gamma);
        CHECK(u[hd_core::RHO] == 2.0);
        CHECK(near(u[hd_core::ENERGY], 7.5 + speed_squared));
        const auto roundtrip = u.to_Primtivies(gamma);
        for (int k = 0; k < hd_core::NVAR; ++k) CHECK(near(roundtrip[k], w[k]));
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            const double normal_v = w[hd_core::VEL(d)];
            CHECK(u[hd_core::MOMDEN(d)] == 2.0 * normal_v);
            const auto flux = hydro::physical_flux(w, d, gamma);
            CHECK(near(flux[hd_core::RHO], 2.0 * normal_v));
            CHECK(near(flux[hd_core::ENERGY], (10.5 + speed_squared) * normal_v));
            for (int transverse = 0; transverse < hd_core::NDIM; ++transverse)
                CHECK(near(flux[hd_core::MOMDEN(transverse)],
                    2.0 * normal_v * w[hd_core::VEL(transverse)] + (transverse == d ? 3.0 : 0.0)));
        }

        // Value operators must not mutate either operand; compound operators must.
        const auto before = u.data;
        const auto result = ((u + u) - (-u)) * 2.0 / 3.0;
        for (int k = 0; k < hd_core::NVAR; ++k) CHECK(near(result[k], 2.0 * u[k]));
        CHECK(u.data == before);
        auto changed = u;
        CHECK(&(changed += u) == &changed);
        changed -= u;
        changed *= 3.0;
        changed /= 3.0;
        for (int k = 0; k < hd_core::NVAR; ++k) CHECK(near(changed[k], u[k]));
        CHECK((2.0 * u).data == (u * 2.0).data);
        expect_throw<std::domain_error>([&] { changed /= 0.0; });
        for (int k = 0; k < hd_core::NVAR; ++k) CHECK(near(changed[k], u[k]));

        // Tiny direct division must not overflow through an intermediate reciprocal.
        hd_core::Conserved tiny{};
        tiny[0] = 1e-310;
        CHECK(near((tiny / 1e-310)[0], 1.0));
        hd_core::Primitives vacuum{};
        vacuum.validate(gamma);
        CHECK(vacuum.to_Conserved(gamma).data == hd_core::Conserved{}.data);
        for (double f : hydro::physical_flux(vacuum, 0, gamma).data) CHECK(f == 0.0);

        auto invalid = w;
        invalid[hd_core::RHO] = -1.0;
        expect_throw<std::invalid_argument>([&] { invalid.validate(gamma); });
        invalid = w;
        invalid[hd_core::PRE] = 0.0;
        expect_throw<std::invalid_argument>([&] { invalid.validate(gamma); });
        invalid = w;
        invalid[hd_core::VEL(0)] = std::numeric_limits<double>::quiet_NaN();
        expect_throw<std::invalid_argument>([&] { invalid.validate(gamma); });
        expect_throw<std::invalid_argument>([&] { w.validate(1.0); });
        auto bad_energy = u;
        bad_energy[hd_core::ENERGY] = speed_squared;
        expect_throw<std::invalid_argument>([&] { bad_energy.validate(gamma); });
        expect_throw<std::out_of_range>([&] { hydro::physical_flux(w, hd_core::NDIM, gamma); });
    });
}
