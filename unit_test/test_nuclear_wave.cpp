#include "test_support.h"
// PressureWave is private to this translation unit. Exercise it directly while
// the public Nuclear Exact solver remains intentionally disabled.
#include "hydro/riemann/exact.cpp"

namespace
{
using namespace hd_microphysics;
constexpr auto nuclear = Equation_of_State::Nuclear;

void relative(double a, double b, double tolerance = 2e-9)
{
    if (!std::isfinite(a) || !std::isfinite(b)
        || std::abs(a - b) > tolerance * std::max(std::abs(a), std::abs(b)))
    {
        std::cerr << "relative mismatch: " << a << " vs " << b
                  << ", error=" << (a-b)/b << '\n';
        throw std::runtime_error("Nuclear wave reference mismatch");
    }
}

hydro::PressureWave wave(double density_ratio, double gamma, double thermal_ratio)
{
    hd_core::Primitives w{};
    w[hd_core::RHO] = density_ratio * nuclear_density_cgs / physcon::units.unit_density;
    double pc, uc;
    nuclear_cold_state(w[hd_core::RHO], gamma, pc, uc);
    w[hd_core::PRE] = pc * (1.0 + thermal_ratio);
    return hydro::PressureWave(w, gamma, nuclear, 1e-9);
}

void ideal_limit()
{
    // gamma_cold == gamma_thermal, below the transition: an independent
    // gamma-law Riemann oracle applies to both shocks and rarefactions.
    const double g = nuclear_gamma_thermal;
    const auto w = wave(0.01, g, 0.7);
    for (double z : {1e-40, 0.01, 0.5, 1.0-1e-9, 1.0, 1.0+1e-9, 1.01, 2.0, 10.0})
    {
        const double p = w.p0*z;
        const auto [f, df] = w.evaluate(p);
        const auto [rho, c] = w.wave_state_at_pressure(p);
        double expected_rho, expected_f, expected_df;
        if (p <= w.p0)
        {
            const double r = p/w.p0;
            expected_rho = w.rho0*std::pow(r, 1.0/g);
            const double exponent = (g-1.0)/(2.0*g);
            expected_f = 2.0*w.sound_speed/(g-1.0)*std::expm1(exponent*std::log(r));
            expected_df = 1.0/(expected_rho*w.sound_speed*std::pow(r,exponent));
        }
        else
        {
            const double b = (g-1.0)/(g+1.0);
            expected_rho = w.rho0*(z+b)/(b*z+1.0);
            const double root = std::sqrt(2.0/((g+1.0)*w.rho0)/(p+b*w.p0));
            expected_f = (p-w.p0)*root;
            expected_df = root*(1.0-0.5*(p-w.p0)/(p+b*w.p0));
        }
        relative(rho, expected_rho);
        relative(c*c, g*p/rho);
        if (p == w.p0) CHECK(f == 0.0); else relative(f, expected_f, 2e-7);
        relative(df, expected_df);
    }
    relative(w.evaluate(0.0).first, -2.0*w.sound_speed/(g-1.0));
    CHECK(w.evaluate(0.0).second == 0.0);
    relative(w.vacuum_limit(), w.evaluate(0.0).first);
}

void cold_transition()
{
    const double g = 4.0/3.0;
    const auto w = wave(2.0, g, 0.0);
    const double cn_high = std::sqrt(nuclear_gamma_stiff*w.p_n/w.rho_n);
    const double cn_low = std::sqrt(g*w.p_n/w.rho_n);
    const double to_transition = 2.0*(cn_high-w.sound_speed)/(nuclear_gamma_stiff-1.0);
    for (double z : {1.1, 1.0, 0.99, 0.3, 1e-12})
    {
        const double p = w.p_n*z;
        const double exponent = z > 1.0 ? nuclear_gamma_stiff : g;
        const double density = w.rho_n*std::pow(z,1.0/exponent);
        const double sound = std::sqrt(exponent*p/density);
        const double reference = z > 1.0
            ? 2.0*(sound-w.sound_speed)/(exponent-1.0)
            : to_transition+2.0*(sound-cn_low)/(g-1.0);
        relative(w.evaluate(p).first, reference);
    }
    relative(w.evaluate(0.0).first, to_transition-2.0*cn_low/(g-1.0));
    const auto at_kink = wave(1.0, g, 0.0);
    relative(at_kink.evaluate(0.5*at_kink.p0).first,
        2.0*at_kink.sound_speed/(g-1.0)*std::expm1((g-1.0)/(2*g)*std::log(0.5)));
}

void hybrid_shocks()
{
    for (double initial_density : {0.5, 0.99, 1.0, 2.0})
    for (double thermal : {0.0, 0.8})
    {
        const auto w = wave(initial_density, 4.0/3.0, thermal);
        for (double compression : {1.0001, 1.02, 1.5, 3.0})
        {
            // Prescribe density, solve the Hugoniot explicitly for pressure.
            // This avoids using the production density root as its own oracle.
            const double rho = w.rho0*compression;
            double pc, uc;
            nuclear_cold_state(rho, w.gamma, pc, uc);
            const double d = (rho-w.rho0)/w.rho0/rho;
            const double a = nuclear_gamma_thermal-1.0;
            const double p = (w.e0-uc/rho+pc/(a*rho)+0.5*w.p0*d)
                /(1.0/(a*rho)-0.5*d);
            CHECK(p > w.p0);
            const auto state = w.wave_state_at_pressure(p);
            relative(state.first, rho);
            const auto [f, df] = w.evaluate(p);
            relative(f, std::sqrt((p-w.p0)*d));
            const double h = 1e-5*(p-w.p0);
            const double finite_difference = (w.evaluate(p+h).first-w.evaluate(p-h).first)/(2*h);
            relative(df, finite_difference, 2e-5);
        }
        for (double dp : {1e-10, 1e-7})
        {
            const double p = w.p0*(1.0+dp);
            const auto [f, df] = w.evaluate(p);
            const double initial_sound = w.wave_state_at_pressure(w.p0, initial_density == 1.0).second;
            relative(f/(p-w.p0), 1.0/w.rho0/initial_sound, 1e-6);
            CHECK(std::isfinite(df) && df > 0.0);
        }
        // Derivative of the integrated rarefaction, including paths that cross
        // the nuclear threshold, checked away from the derivative discontinuity.
        for (double z : {0.1, 0.7})
        {
            const double p = w.p0*z, h = 1e-5*p;
            const auto [f, df] = w.evaluate(p);
            CHECK(f < 0.0);
            relative(df, (w.evaluate(p+h).first-w.evaluate(p-h).first)/(2*h), 2e-7);
        }
    }
}
}

int main()
{
    return run_test("Nuclear PressureWave", [] {
        for (bool scaled : {false, true})
        {
            if (scaled) physcon::units.init_unitsystem_RHOLV(1e12,1e6,1e9,1e10);
            else physcon::units.init_unitsystem_MLT(1,1,1,1);
            ideal_limit();
            cold_transition();
            hybrid_shocks();
            expect_throw<std::invalid_argument>([] { wave(0.5,4.0/3.0,-0.1); });
        }
        physcon::units.init_unitsystem_MLT(1,1,1,1);
    });
}
