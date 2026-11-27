#include "test_support.h"
// Include the implementation to test the private PressureWave against analytic curves.
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
        hd_core::Primitives upstream{};
        upstream[hd_core::RHO] = w.rho0;
        upstream[hd_core::PRE] = w.p0;
        const int normal = hd_core::NDIM - 1;
        const double u0 = 0.3*w.sound_speed;
        upstream[hd_core::VEL(normal)] = u0;
        auto check_outer_speed = [&](double p, double relative_speed) {
            relative(hydro::outer_wave_speed(upstream, w, p, hydro::Side::Left, normal),
                     u0 - relative_speed);
            relative(hydro::outer_wave_speed(upstream, w, p, hydro::Side::Right, normal),
                     u0 + relative_speed);
        };
        // The rarefaction head travels in the original gas, even across the kink or into vacuum.
        for (double p : {0.0, 0.1*w.p0, w.p0})
            check_outer_speed(p, w.sound_speed);
        if (initial_density == 1.0)
        {
            const double c_high = std::sqrt((nuclear_gamma_stiff*w.p_n
                + nuclear_gamma_thermal*w.p_th0) / w.rho_n);
            check_outer_speed(std::nextafter(w.p0, std::numeric_limits<double>::infinity()), c_high);
        }
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
            check_outer_speed(p, std::sqrt((p-w.p0)/d)/w.rho0);
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

void pressure_brackets()
{
    // Below nuclear density, equal cold/thermal exponents give an
    // independent gamma-law oracle for the complete star-pressure solve.
    const double g = nuclear_gamma_thermal;
    const auto left = wave(0.001, g, 0.7);
    const auto right = wave(0.002, g, 0.7);
    hd_core::Primitives L{}, R{};
    L[hd_core::RHO] = left.rho0; L[hd_core::PRE] = left.p0;
    R[hd_core::RHO] = right.rho0; R[hd_core::PRE] = right.p0;
    auto reference = [&](const hydro::PressureWave& w, double p) {
        if (p <= w.p0)
            return 2.0*w.sound_speed/(g-1.0)
                * std::expm1((g-1.0)/(2*g)*std::log(p/w.p0));
        return (p-w.p0)*std::sqrt(2.0/((g+1.0)*w.rho0)
                /(p+(g-1.0)/(g+1.0)*w.p0));
    };
    // Positive lower bounds, bracket expansion, and the vacuum-limit fallback.
    for (double ratio : {0.25, 1.5, 8.0, 1e-4})
    {
        const double p = ratio*left.p0;
        const double fL = reference(left, p), fR = reference(right, p);
        R[hd_core::VEL(0)] = -fL-fR;
        for (auto method : {hydro::ExactPressureMethod::SafeguardedNewton,
                            hydro::ExactPressureMethod::Bisection})
        {
            const auto star = hydro::solve_star(L, R, left, right, 1e-9,
                                               200, 0, method, 0.0);
            CHECK(star.has_value());
            relative(star->first, p, 2e-8);
            relative(star->second, -fL, 2e-8);
        }
    }
    R[hd_core::VEL(0)] = 2.1*(left.sound_speed+right.sound_speed)/(g-1.0);
    for (auto method : {hydro::ExactPressureMethod::SafeguardedNewton,
                        hydro::ExactPressureMethod::Bisection})
        CHECK(!hydro::solve_star(L, R, left, right, 1e-9, 200, 0, method, 0.0));
}

void near_transition_integrals()
{
    for (double density : {1.0+1e-12, 1.0001, 1.01, 1.2, 2.0})
    for (double thermal : {0.01, 0.8})
    {
        const auto w = wave(density, 4.0/3.0, thermal);
        const double pt = w.p_n + w.p_th0*std::pow(w.rho_n/w.rho0, nuclear_gamma_thermal);
        const double ft = w.evaluate(pt).first;
        for (double offset : {-1e-10, -1e-12, -1e-14, 1e-14, 1e-12, 1e-10})
        {
            const double p = pt*(1.0+offset);
            try
            {
                const auto [f, df] = w.evaluate(p);
                CHECK(std::isfinite(f) && std::isfinite(df));
                const double c = std::sqrt(((offset > 0 ? nuclear_gamma_stiff : w.gamma)*w.p_n
                                  + nuclear_gamma_thermal*(pt-w.p_n))/w.rho_n);
                relative(df, 1.0/(w.rho_n*c), 2e-9);
                CHECK(std::abs(f-ft) <= 2*std::abs(p-pt)/(w.rho_n*c)
                      + 1e-9*std::abs(ft));
            }
            catch (...)
            {
                std::cerr << "transition density=" << density << " thermal=" << thermal
                          << " offset=" << offset << '\n';
                throw;
            }
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
            pressure_brackets();
            near_transition_integrals();
            expect_throw<std::invalid_argument>([] { wave(0.5,4.0/3.0,-0.1); });
        }
        physcon::units.init_unitsystem_MLT(1,1,1,1);
    });
}
