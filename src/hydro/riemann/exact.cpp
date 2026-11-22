/**
 * @file exact.cpp
 * @brief Hydrodynamic spatial operators.
 */
#include "exact.h"
#include "utils/performance.h"
#include "core/cell.h"
#include "microphysics/eos/eos.h"
#include "microphysics/composition/species.h"
#include "physcon/constants.h"
#include "flux.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>
namespace hydro
{
    namespace
    {
        /** @brief Left or right side of a Riemann problem. */
        enum class Side
        {
            Left, ///< Left side of the interface.
            Right ///< Right side of the interface.
        };
        // Aggregate one solve on exit, including failures. No clock reads or allocation.
        struct RootSolveRecord
        {
            std::uint64_t iterations = 0, evaluations = 0;
            std::uint64_t initial_guesses = 0, newton_steps = 0;
            std::uint64_t bisection_steps = 0, bracket_expansions = 0;
            bool failed = true, vacuum = false;
            ~RootSolveRecord()
            {
                if (!performance::enabled()) return;
                auto& r = performance::exact_roots;
                ++r.solves;
                r.iterations += iterations;
                r.evaluations += evaluations;
                r.initial_guesses += initial_guesses;
                r.newton_steps += newton_steps;
                r.bisection_steps += bisection_steps;
                r.bracket_expansions += bracket_expansions;
                r.failures += failed;
                r.vacuum += vacuum;
                r.max_iterations = std::max(r.max_iterations, iterations);
                const auto bin = iterations == 0 ? 0 : iterations <= 4 ? 1 : iterations <= 8 ? 2 :
                    iterations <= 16 ? 3 : iterations <= 32 ? 4 : iterations <= 64 ? 5 : 6;
                ++r.histogram[bin];
            }
        };
        struct WavePoint
        {
            double value, derivative; // f(p), df/dp; u*_L=u_L-f_L, u*_R=u_R+f_R.
            double density, sound_speed;
        };

        // One side's ORIGINAL state, unchanged throughout the star-pressure solve.
        // Composition/Saha are frozen. Only the radiation branch needs mu, obtained
        // once from the incoming face state. No composition scan occurs in a wave
        // iteration. Constant specific ionization energy cancels in the Hugoniot.
        struct PressureWave
        {
            double rho0, pressure, gamma, sound_speed;
            double shock_a, shock_b;
            double rarefaction_exponent, rarefaction_coefficient;
            bool radiative;
            double gas_fraction = 1.0, radiation_fraction = 0.0;
            double energy_scale, velocity_scale, e0, inner_tolerance;
            mutable std::optional<double> vacuum_value;

            // Solve b*t + d*t^4 = 1 with 0<t<=1, b,d<=1 and b+d>=1.
            static double scaled_temperature(double b, double d)
            {
                double lo = 0.0, hi = 1.0, t = 1.0;
                for (int n = 0; n < 80; ++n)
                {
                    const double t3 = t*t*t;
                    const double f = b*t + d*t*t3 - 1.0;
                    if (std::abs(f) <= 8.0 * std::numeric_limits<double>::epsilon())
                        return t;
                    if (f > 0.0) hi = t; else lo = t;
                    const double next = t - f / (b + 4.0*d*t3);
                    t = std::isfinite(next) && next > lo && next < hi
                        ? next : lo + 0.5*(hi-lo);
                }
                throw std::runtime_error("ExactSolver: radiation temperature did not converge");
            }

            PressureWave(const hd_core::Primitives& W, double g,
                         hd_microphysics::Equation_of_State model, double tolerance)
                : rho0(W[hd_core::RHO]), pressure(W[hd_core::PRE]), gamma(g),
                  sound_speed(std::sqrt(g * pressure / rho0)),
                  shock_a(2 / ((g + 1) * rho0)),
                  shock_b((g - 1) / (g + 1) * pressure),
                  rarefaction_exponent((g - 1) / (2 * g)),
                  rarefaction_coefficient(2 * sound_speed / (g - 1)),
                  radiative(model == hd_microphysics::Equation_of_State::Ionized_Radiative),
                  energy_scale(pressure/rho0), velocity_scale(std::sqrt(energy_scale)),
                  e0(energy_scale/(g-1)),
                  inner_tolerance(std::clamp(0.01*tolerance, 5e-14, 1e-10))
            {
                if (!std::isfinite(sound_speed) || sound_speed <= 0.0)
                    throw std::invalid_argument("ExactSolver: invalid initial sound speed");
                if (!radiative) return;
                hd_core::Cell context;
                context.U[hd_core::RHO] = rho0;
                for (int k = 0; k < hd_core::NCOMP; ++k)
                    context.U[hd_core::COMP(k)] = rho0 * W[hd_core::COMP(k)];
                hd_microphysics::validate_cell_composition(context);
                hd_microphysics::evaluate_cell_density(context);
                const double R = physcon::R_g / context.mean_molecular_weight;
                if (!std::isfinite(R) || R <= 0.0 || !std::isfinite(physcon::a_Rad)
                    || physcon::a_Rad <= 0.0)
                    throw std::invalid_argument("ExactSolver: invalid radiation EoS parameters");

                // Scale by the smaller gas-only / radiation-only temperature.
                // Logarithms avoid forming T^4 or rho*R at an extreme scale.
                const double log_b = std::log(rho0) + std::log(R);
                const double log_d = std::log(physcon::a_Rad/3.0);
                const double log_p = std::log(pressure);
                const bool gas_limited = log_p-log_b <= 0.25*(log_p-log_d);
                const double log_T = gas_limited ? log_p-log_b : 0.25*(log_p-log_d);
                // The limiting coefficient is exactly one, preserving the
                // temperature bracket even if large logarithms lose a few ulps.
                const double b = gas_limited ? 1.0 : std::exp(std::min(0.0, log_b+log_T-log_p));
                const double d = gas_limited ? std::exp(std::min(0.0, log_d+4.0*log_T-log_p)) : 1.0;
                const double t = scaled_temperature(b, d);
                gas_fraction = b*t;
                radiation_fraction = d*t*t*t*t;
                e0 = energy_scale * (gas_fraction/(g-1) + 3.0*radiation_fraction);
                sound_speed = velocity_scale * std::sqrt(thermo(1.0, 1.0).cs2);
                if (!std::isfinite(e0) || !std::isfinite(sound_speed) || sound_speed <= 0.0)
                    throw std::runtime_error("ExactSolver: initial radiation state exceeds numerical range");
            }

            // Dimensionless rho/rho0, P/P0, e/(P0/rho0), c^2/(P0/rho0).
            // e contains thermal+radiation energy only; e_ion(X) cancels exactly.
            struct Thermo { double e, cs2, e_rho, e_p, beta; };
            Thermo thermo(double r, double p) const
            {
                if (!std::isfinite(r) || !std::isfinite(p) || r <= 0.0 || p <= 0.0)
                    throw std::runtime_error("ExactSolver: wave left the positive EoS domain");
                const double log_p = std::log(p);
                const double log_b = std::log(gas_fraction) + std::log(r);
                const double log_d = std::log(radiation_fraction);
                const bool gas_limited = log_p-log_b <= 0.25*(log_p-log_d);
                const double log_t = gas_limited ? log_p-log_b : 0.25*(log_p-log_d);
                const double b = gas_limited ? 1.0 : std::exp(std::min(0.0, log_b+log_t-log_p));
                const double d = gas_limited ? std::exp(std::min(0.0, log_d+4.0*log_t-log_p)) : 1.0;
                const double t = scaled_temperature(b, d);
                // Compute pressures as fractions of p, not by raising T to 4.
                const double pg = p*b*t;
                const double pr = p*d*t*t*t*t;
                const double A = pg + 4.0*pr;
                const double B = pg/(gamma-1.0) + 12.0*pr;
                const double e = (pg/(gamma-1.0) + 3.0*pr)/r;
                const double cs2 = (pg + A*(A/B))/r;
                // T*e_T = B/r and T*P_T = A.
                const double e_p = (B/A)/r;
                const double e_rho = -(3.0*pr/r + (B/A)*(pg/r))/r;
                if (!std::isfinite(e) || !std::isfinite(cs2) || cs2 <= 0.0)
                    throw std::runtime_error("ExactSolver: nonfinite radiation wave state");
                return {e, cs2, e_rho, e_p, pg/(pg+pr)};
            }

            // Bound Gamma1 along the remaining expansion. For gamma<4/3,
            // beta=Pg/P decreases toward zero; for gamma>4/3 it increases to 1.
            // Gamma1(beta) is convex, but can have an interior minimum: simply
            // using min(gamma,4/3) is not a valid bound for every gas gamma.
            double minimum_gamma1(double beta) const
            {
                double lo = gamma < 4.0/3.0 ? 0.0 : beta;
                double hi = gamma < 4.0/3.0 ? beta : 1.0;
                for (int n = 0; n < 64; ++n)
                {
                    const double b = lo+0.5*(hi-lo), A = 4.0-3.0*b;
                    const double B = b/(gamma-1.0)+12.0*(1.0-b);
                    const double derivative = 1.0-6.0*A/B
                        - (A/B)*(A/B)*(1.0/(gamma-1.0)-12.0);
                    if (derivative > 0.0) hi = b; else lo = b;
                }
                const double b = lo+0.5*(hi-lo), A = 4.0-3.0*b;
                const double B = b/(gamma-1.0)+12.0*(1.0-b);
                return b+A*A/B-64.0*std::numeric_limits<double>::epsilon();
            }

            using Integral = std::array<double, 2>; // {log(rho/rho0), f/velocity_scale}
            Integral rarefaction_rhs(double q, const Integral& y) const
            {
                const double r = std::exp(y[0]), p = std::exp(q);
                const auto s = thermo(r, p);
                return {p/r/s.cs2, p/r/std::sqrt(s.cs2)};
            }

            Integral rk4(double q, const Integral& y, double h) const
            {
                auto shifted = [&](const Integral& k, double factor) {
                    return Integral{y[0]+factor*h*k[0], y[1]+factor*h*k[1]};
                };
                const auto a = rarefaction_rhs(q, y);
                const auto b = rarefaction_rhs(q+0.5*h, shifted(a, 0.5));
                const auto c = rarefaction_rhs(q+0.5*h, shifted(b, 0.5));
                const auto d = rarefaction_rhs(q+h, shifted(c, 1.0));
                return {y[0]+h/6.0*(a[0]+2*b[0]+2*c[0]+d[0]),
                        y[1]+h/6.0*(a[1]+2*b[1]+2*c[1]+d[1])};
            }

            // Integrate in log pressure: neither rho nor P can become negative.
            // Step doubling controls both log density and the velocity integral.
            Integral integrate(double from, double to, Integral y) const
            {
                if (to > from || to < -640.0)
                    throw std::runtime_error("ExactSolver: rarefaction exceeds numerical range");
                double q = from, h = -0.25;
                for (int n = 0; q > to; ++n)
                {
                    if (n >= 50000)
                        throw std::runtime_error("ExactSolver: rarefaction integration did not converge");
                    h = std::max(h, to-q);
                    if (q+h == q)
                        throw std::runtime_error("ExactSolver: rarefaction step underflow");
                    const auto coarse = rk4(q, y, h);
                    const auto half = rk4(q, y, 0.5*h);
                    const auto fine = rk4(q+0.5*h, half, 0.5*h);
                    const double error = std::max(std::abs(fine[0]-coarse[0]),
                        std::abs(fine[1]-coarse[1])/std::max(1.0, std::abs(fine[1]))) / 15.0;
                    if (!std::isfinite(error))
                        throw std::runtime_error("ExactSolver: nonfinite rarefaction integration error");
                    if (error <= inner_tolerance)
                    {
                        for (int j = 0; j < 2; ++j)
                            y[j] = fine[j] + (fine[j]-coarse[j])/15.0;
                        q = h == to-q ? to : q+h;
                    }
                    const double factor = error > 0.0
                        ? std::clamp(0.9*std::pow(inner_tolerance/error, 0.2), 0.2, 2.0)
                        : 2.0;
                    h = -std::min(4.0, std::abs(h)*factor);
                }
                return y;
            }

            double vacuum_limit() const
            {
                if (!radiative) return -rarefaction_coefficient;
                if (vacuum_value) return *vacuum_value;
                Integral y{0.0, 0.0};
                double q = 0.0;
                for (int n = 0; n < 160; ++n)
                {
                    y = integrate(q, q-4.0, y);
                    q -= 4.0;
                    const auto s = thermo(std::exp(y[0]), std::exp(q));
                    const double c = std::sqrt(s.cs2), gmin = minimum_gamma1(s.beta);
                    // Gamma1>=gmin>1 bounds the remaining integral by
                    // 2*c/(gmin-1). Never pass P=0 to the thermal EoS.
                    if (gmin > 1.0 && 2.0*c/(gmin-1.0)
                        <= inner_tolerance*std::max(1.0, sound_speed/velocity_scale))
                    {
                        vacuum_value = velocity_scale*y[1];
                        return *vacuum_value;
                    }
                }
                throw std::runtime_error("ExactSolver: vacuum limit did not converge within numerical range");
            }

            WavePoint radiation_point(double p) const
            {
                if (p == 0.0) return {vacuum_limit(), 0.0, 0.0, 0.0};
                const double z = p/pressure, dz = (p-pressure)/pressure;
                if (p == pressure)
                    return {0.0, 1.0/(rho0*sound_speed), rho0, sound_speed};
                // Shock and isentropic curves agree through second order at the
                // original state. Use that limit to avoid cancellation in H and
                // log(P/P0); retain the first-order change in the shock speed.
                if (std::abs(dz) < 1e-7)
                {
                    const double beta = gas_fraction/(gas_fraction+radiation_fraction);
                    const double A = 4.0-3.0*beta;
                    const double B = beta/(gamma-1.0)+12.0*(1.0-beta);
                    const double G1 = beta+A*A/B;
                    const double dG1_db = 1.0-6.0*A/B
                        - (A/B)*(A/B)*(1.0/(gamma-1.0)-12.0);
                    const double db_dlnrho = -beta*(1.0-beta)*(3.0*A/B-1.0);
                    const double fundamental = 1.0+0.5*(G1-1.0+dG1_db*db_dlnrho/G1);
                    const double x = (p-pressure)/(rho0*sound_speed*sound_speed);
                    const double r = 1.0+x-(fundamental-1.0)*x*x;
                    const double c = velocity_scale*std::sqrt(thermo(r, z).cs2);
                    return {(p-pressure)/(rho0*sound_speed)*(1.0-0.5*fundamental*x),
                            (1.0-fundamental*x)/(rho0*sound_speed), rho0*r, c};
                }
                if (p < pressure)
                {
                    const double q = std::log(p)-std::log(pressure);
                    const auto y = integrate(0.0, q, {0.0, 0.0});
                    const double r = std::exp(y[0]);
                    const double c = velocity_scale*std::sqrt(thermo(r, z).cs2);
                    return {velocity_scale*y[1], 1.0/(rho0*r*c), rho0*r, c};
                }

                // H(r,z)=e(r,z)-e0-(z+1)/2*(1-1/r)=0, r>=1.
                const double initial_e = e0/energy_scale;
                auto H = [&](double r) {
                    return thermo(r, z).e - initial_e - 0.5*(z+1.0)*((r-1.0)/r);
                };
                double lo = 1.0, hi = 2.0;
                for (int n = 0; H(hi) > 0.0; ++n)
                {
                    if (n >= 100 || !std::isfinite(hi*2.0))
                        throw std::runtime_error("ExactSolver: cannot bracket shock density");
                    hi *= 2.0;
                }
                // Resolve the compression itself, also for weak shocks.
                for (int n = 0; n < 160; ++n)
                {
                    const double r = lo+0.5*(hi-lo);
                    if (r == lo || r == hi || hi-lo <= inner_tolerance*(r-1.0)) break;
                    if (H(r) > 0.0) lo = r; else hi = r;
                    if (n == 159)
                        throw std::runtime_error("ExactSolver: shock density did not converge");
                }
                const double r = lo+0.5*(hi-lo), D = (r-1.0)/r;
                const auto s = thermo(r, z);
                const double f = std::sqrt(dz*D);
                const double dr_dp = (0.5*D-s.e_p)/(s.e_rho-0.5*(z+1.0)/r/r);
                const double df_dp = (D+dz*dr_dp/r/r)/(2.0*f);
                if (!std::isfinite(df_dp) || df_dp <= 0.0)
                    throw std::runtime_error("ExactSolver: invalid shock pressure derivative");
                return {velocity_scale*f, velocity_scale/pressure*df_dp,
                        rho0*r, velocity_scale*std::sqrt(s.cs2)};
            }

            // Return {wave value, pressure derivative}, sharing sqrt/pow work.
            std::pair<double, double> evaluate(double trial_pressure) const
            {
                if (!std::isfinite(trial_pressure) || trial_pressure < 0.0)
                    throw std::invalid_argument("ExactSolver: invalid trial pressure");
                if (radiative)
                {
                    const auto s = radiation_point(trial_pressure);
                    return {s.value, s.derivative};
                }
                if (trial_pressure > pressure)
                {
                    const double root = std::sqrt(shock_a / (trial_pressure + shock_b));
                    return {
                        (trial_pressure - pressure) * root,
                        root * (1 - (trial_pressure - pressure) / (2 * (trial_pressure + shock_b)))
                    };
                }
                const double power = std::pow(trial_pressure / pressure, rarefaction_exponent);
                return {
                    rarefaction_coefficient * (power - 1),
                    trial_pressure > 0
                        ? rarefaction_coefficient * rarefaction_exponent * power / trial_pressure
                        : 0
                };
            }
        };

        // Scalar adapter: both root finders now use the same wave relation.
        double wave_condition(double p, const PressureWave& wave)
        {
            return wave.evaluate(p).first;
        }

        double outer_wave_speed(const hd_core::Primitives& W, const PressureWave& wave,
                                double p_star, Side side, int direction)
        {
            double speed = wave.sound_speed;
            if (p_star > wave.pressure)
            {
                if (wave.radiative)
                {
                    const double dp = p_star-wave.pressure;
                    // j=(P*-P0)/f, S-u0=sign*j/rho0; avoids subtracting
                    // two almost equal densities for a very weak shock.
                    speed = (dp/wave.rho0)/wave.radiation_point(p_star).value;
                }
                else
                {
                    const double g = wave.gamma;
                    speed *= std::sqrt((g+1.0)/(2.0*g)*(p_star/wave.pressure)
                                     + (g-1.0)/(2.0*g));
                }
            }
            return W[hd_core::VEL(direction)] + (side == Side::Left ? -speed : speed);
        }

        hd_core::Primitives sample_rarefraction(
            const hd_core::Primitives& K, const PressureWave& wave, double v_star,
            double p_star, Side side, int direction)
        {
            const double sign = side == Side::Left ? -1.0 : 1.0;
            const double u0 = K[hd_core::VEL(direction)];
            if (sign*(u0 + sign*wave.sound_speed) <= 0.0) return K;
            hd_core::Primitives face = K; // X and tangential velocity are advected.
            if (!wave.radiative)
            {
                const double g = wave.gamma;
                const double c_star = wave.sound_speed
                    * std::pow(p_star/wave.pressure, (g-1.0)/(2.0*g));
                if (sign*(v_star + sign*c_star) >= 0.0)
                {
                    face[hd_core::RHO] = wave.rho0*std::pow(p_star/wave.pressure, 1.0/g);
                    face[hd_core::VEL(direction)] = v_star;
                    face[hd_core::PRE] = p_star;
                }
                else
                {
                    const double u = ((g-1.0)*u0-sign*2.0*wave.sound_speed)/(g+1.0);
                    const double c = -sign*u;
                    face[hd_core::RHO] = wave.rho0*std::pow(c/wave.sound_speed, 2.0/(g-1.0));
                    face[hd_core::VEL(direction)] = u;
                    face[hd_core::PRE] = wave.pressure*std::pow(c/wave.sound_speed, 2.0*g/(g-1.0));
                }
                return face;
            }

            const auto star = wave.radiation_point(p_star);
            if (sign*(v_star + sign*star.sound_speed) >= 0.0)
            {
                face[hd_core::RHO] = star.density;
                face[hd_core::VEL(direction)] = v_star;
                face[hd_core::PRE] = p_star;
                return face;
            }
            // The face lies inside the fan. Solve u(P)+sign*c(P)=x/t=0.
            // Use log pressure so near-vacuum fans do not need linear-P bisection.
            auto characteristic = [&](double q) {
                const auto s = wave.radiation_point(wave.pressure*std::exp(q));
                return sign*u0 + s.value + s.sound_speed;
            };
            double lo = p_star > 0.0 ? std::log(p_star)-std::log(wave.pressure) : -1.0;
            double hi = 0.0;
            if (p_star == 0.0)
            {
                while (characteristic(lo) > 0.0)
                {
                    if (lo <= -640.0)
                        throw std::runtime_error("ExactSolver: cannot bracket rarefaction fan sample");
                    lo = std::max(-640.0, 2.0*lo);
                }
            }
            for (int n = 0; n < 100; ++n)
            {
                const double q = lo+0.5*(hi-lo);
                if (q == lo || q == hi || hi-lo <= wave.inner_tolerance) break;
                if (characteristic(q) > 0.0) hi = q; else lo = q;
                if (n == 99)
                    throw std::runtime_error("ExactSolver: rarefaction fan sample did not converge");
            }
            const double p = wave.pressure*std::exp(lo+0.5*(hi-lo));
            const auto sample = wave.radiation_point(p);
            face[hd_core::RHO] = sample.density;
            face[hd_core::VEL(direction)] = u0 + sign*sample.value;
            face[hd_core::PRE] = p;
            return face;
        }

        hd_core::Primitives sample_shock(
            const hd_core::Primitives& K, const PressureWave& wave, double v_star,
            double p_star, Side side, int direction)
        {
            const double speed = outer_wave_speed(K, wave, p_star, side, direction);
            if ((side == Side::Left && speed > 0.0) || (side == Side::Right && speed < 0.0))
                return K;
            hd_core::Primitives face = K;
            const double b = (wave.gamma-1.0)/(wave.gamma+1.0);
            const double ratio = p_star/wave.pressure;
            face[hd_core::RHO] = wave.radiative ? wave.radiation_point(p_star).density
                : wave.rho0*(ratio+b)/(b*ratio+1.0);
            face[hd_core::VEL(direction)] = v_star;
            face[hd_core::PRE] = p_star;
            return face;
        }

        InterfaceState sample_interface(
            const hd_core::Primitives& L, const hd_core::Primitives& R,
            const PressureWave& left, const PressureWave& right,
            double v_star, double p_star, int direction)
        {
            const bool from_left = v_star > 0.0;
            const auto& K = from_left ? L : R;
            const auto& wave = from_left ? left : right;
            const Side side = from_left ? Side::Left : Side::Right;
            return {p_star > wave.pressure
                ? sample_shock(K, wave, v_star, p_star, side, direction)
                : sample_rarefraction(K, wave, v_star, p_star, side, direction), wave.gamma};
        }

        // The optimized path uses cached wave coefficients and analytic derivatives.
        std::optional<std::pair<double, double>> solve_star_newton(
            const hd_core::Primitives& W_L, const hd_core::Primitives& W_R,
            const PressureWave& left, const PressureWave& right, double tolerance,
            int max_iterations, int direction, double pressure_hint)
        {
            RootSolveRecord record;
            const double rho_L = W_L[hd_core::RHO], rho_R = W_R[hd_core::RHO];
            const double pre_L = W_L[hd_core::PRE], pre_R = W_R[hd_core::PRE];
            const double v_nL = W_L[hd_core::VEL(direction)], v_nR = W_R[hd_core::VEL(direction)];
            if (pre_L == pre_R && v_nL == v_nR)
            {
                record.failed = false;
                return std::pair<double, double>{pre_L, v_nL};
            }
            auto residual = [&](double p) {
                ++record.evaluations;
                const auto l = left.evaluate(p), r = right.evaluate(p);
                return std::pair<double, double>{l.first + r.first + v_nR - v_nL, l.second + r.second};
            };
            if (residual(0).first >= 0)
            {
                record.failed = false;
                record.vacuum = true;
                return std::nullopt;
            }
            double low = 0, high = std::max(pre_L, pre_R);
            while (residual(high).first < 0)
            {
                high *= 2;
                ++record.bracket_expansions;
                if (!std::isfinite(high))
                    throw std::runtime_error("ExactSolver: pressure bracket overflow");
            }
            double p = (low + high) / 2;
            // A cached pressure is only an initial iterate: bracket construction,
            // exact residuals, convergence certification and fallback are unchanged.
            if (std::isfinite(pressure_hint) && pressure_hint > low && pressure_hint <= high)
                p = pressure_hint;
            else
            {
                // Acoustic (PVRS) estimate is inexpensive and accurate for smooth flow.
                // It is only a starting guess; all iterates still use the exact equation.
                const double guess = 0.5 * (pre_L + pre_R) - 0.125 * (v_nR - v_nL) *
                    (rho_L + rho_R) * (left.sound_speed + right.sound_speed);
                if (std::isfinite(guess) && guess > low && guess <= high) p = guess;
            }
            bool current_is_newton = false;
            bool converged = false;
            for (int iteration = 0; iteration < max_iterations; ++iteration)
            {
                ++record.iterations;
                if (iteration == 0) ++record.initial_guesses;
                else if (current_is_newton) ++record.newton_steps;
                else ++record.bisection_steps;
                const auto [f, derivative] = residual(p);
                if (!std::isfinite(f))
                    throw std::runtime_error("ExactSolver: nonfinite pressure residual");
                if (f > 0) high = p;
                else low = p;
                if (high - low <= high * tolerance)
                {
                    converged = true;
                    break;
                }
                const double correction = f / derivative;
                if (std::isfinite(correction) && std::abs(correction) <= p * tolerance)
                {
                    // Certify a narrow sign-changing bracket: a small Newton step alone
                    // is insufficient near vacuum or when cancellation masks the residual.
                    const double lo = std::max(low, p * (1 - 0.25 * tolerance));
                    const double hi = std::min(high, p * (1 + 0.25 * tolerance));
                    if (lo > 0 && residual(lo).first <= 0 && residual(hi).first >= 0 &&
                        hi - lo <= hi * tolerance)
                    {
                        low = lo;
                        high = hi;
                        converged = true;
                        break;
                    }
                }
                const double candidate = p - correction;
                // Force a bisection periodically to bound stagnation for difficult states.
                if (iteration % 8 != 7 && std::isfinite(candidate) &&
                    candidate > low && candidate < high)
                {
                    p = candidate;
                    current_is_newton = true;
                }
                else
                {
                    p = (low + high) / 2;
                    current_is_newton = false;
                }
            }
            if (!converged)
                throw std::runtime_error("ExactSolver: star pressure did not converge");
            const double pre_star = (high + low) / 2;
            const double f_left = left.evaluate(pre_star).first;
            const double f_right = right.evaluate(pre_star).first;
            const double v_nstar = (v_nL - f_left + v_nR + f_right) / 2;
            record.failed = false;
            return std::pair<double, double>{pre_star, v_nstar};
        }

        // Reference bisection uses exactly the same wave physics as Newton.
        std::optional<std::pair<double, double>> solve_star_bisection(
            const hd_core::Primitives& W_L, const hd_core::Primitives& W_R,
            const PressureWave& left, const PressureWave& right, double tolerance,
            int max_iterations, int direction)
        {
            RootSolveRecord record;
            const double pre_L = W_L[hd_core::PRE], pre_R = W_R[hd_core::PRE];
            const double v_nL = W_L[hd_core::VEL(direction)], v_nR = W_R[hd_core::VEL(direction)];
            if (pre_L == pre_R && v_nL == v_nR)
            {
                record.failed = false;
                return std::pair<double, double>{pre_L, v_nL};
            }
            auto residual = [&](double p) {
                ++record.evaluations;
                return wave_condition(p, left) +
                    wave_condition(p, right) + v_nR - v_nL;
            };
            if (residual(0) >= 0)
            {
                record.failed = false;
                record.vacuum = true;
                return std::nullopt;
            }
            double low = 0, high = std::max(pre_L, pre_R);
            while (residual(high) < 0)
            {
                high *= 2;
                ++record.bracket_expansions;
                if (!std::isfinite(high))
                    throw std::runtime_error("ExactSolver: pressure bracket overflow");
            }
            double p = (low + high) / 2;
            bool converged = false;
            for (int iteration = 0; iteration < max_iterations; ++iteration)
            {
                ++record.iterations;
                if (iteration == 0) ++record.initial_guesses;
                else ++record.bisection_steps;
                const double f = residual(p);
                if (!std::isfinite(f))
                    throw std::runtime_error("ExactSolver: nonfinite pressure residual");
                if (f > 0) high = p;
                else low = p;
                if (high - low <= high * tolerance)
                {
                    converged = true;
                    break;
                }
                p = (low + high) / 2;
            }
            if (!converged)
                throw std::runtime_error("ExactSolver: star pressure did not converge");
            const double pre_star = (high + low) / 2;
            const double f_left = wave_condition(pre_star, left);
            const double f_right = wave_condition(pre_star, right);
            const double v_nstar = (v_nL - f_left + v_nR + f_right) / 2;
            record.failed = false;
            return std::pair<double, double>{pre_star, v_nstar};
        }

        std::optional<std::pair<double, double>> solve_star(
            const hd_core::Primitives& W_L, const hd_core::Primitives& W_R,
            const PressureWave& left, const PressureWave& right, double tolerance,
            int max_iterations, int direction, ExactPressureMethod method, double pressure_hint)
        {
            if (method == ExactPressureMethod::Bisection)
                return solve_star_bisection(W_L, W_R, left, right,
                                            tolerance, max_iterations, direction);
            return solve_star_newton(W_L, W_R, left, right,
                                     tolerance, max_iterations, direction, pressure_hint);
        }

        InterfaceState sample_vacuum(
            const hd_core::Primitives& L, const hd_core::Primitives& R,
            const PressureWave& left, const PressureWave& right, int direction)
        {
            const double edge_left = L[hd_core::VEL(direction)] - left.vacuum_limit();
            const double edge_right = R[hd_core::VEL(direction)] + right.vacuum_limit();
            if (edge_left <= 0.0 && edge_right >= 0.0)
                return {hd_core::Primitives{}, left.gamma};
            if (edge_left > 0.0)
                return {sample_rarefraction(L, left, edge_left, 0.0, Side::Left, direction), left.gamma};
            return {sample_rarefraction(R, right, edge_right, 0.0, Side::Right, direction), right.gamma};
        }

        InterfaceState gas_to_vacuum(const hd_core::Primitives& gas,
                                    const PressureWave& wave, Side side, int direction)
        {
            const double sign = side == Side::Left ? -1.0 : 1.0;
            const double edge = gas[hd_core::VEL(direction)] + sign*wave.vacuum_limit();
            if (sign*edge >= 0.0) return {hd_core::Primitives{}, wave.gamma};
            return {sample_rarefraction(gas, wave, edge, 0.0, side, direction), wave.gamma};
        }

    }
    ExactSolver::ExactSolver(double tolerance, int max_iterations, ExactPressureMethod method)
    : tolerance_(tolerance),
      max_iterations_(max_iterations), method_(method)
    {
        if (method_ != ExactPressureMethod::SafeguardedNewton && method_ != ExactPressureMethod::Bisection)
            throw std::invalid_argument("ExactSolver: invalid pressure method");
        if (!std::isfinite(tolerance_) || tolerance_ <= 0.0)
        {
            throw std::invalid_argument(
                "ExactSolver: tolerance must be finite and positive"
            );
        }

        if (max_iterations_ <= 0)
        {
            throw std::invalid_argument(
                "ExactSolver: max_iterations must be positive"
            );
        }
    }
    hd_core::Conserved ExactSolver::flux(
        const hd_core::Primitives& WL, const hd_core::Primitives& WR,
        int direction, double gamma_L, double gamma_R, double& face_speed, double& face_velocity) const
    {
        return flux_impl(WL, WR, direction, gamma_L, gamma_R, nullptr, face_speed, face_velocity);
    }

    hd_core::Conserved ExactSolver::flux_with_pressure_hint(
        const hd_core::Primitives& WL, const hd_core::Primitives& WR,
        int direction, double gamma_L, double gamma_R, double& pressure_hint, double& face_speed, double& face_velocity) const
    {
        return flux_impl(WL, WR, direction, gamma_L, gamma_R, &pressure_hint, face_speed, face_velocity);
    }

        hd_core::Conserved ExactSolver::flux_impl(
            const hd_core::Primitives& WL,
            const hd_core::Primitives& WR,
            int direction,
            double gamma_L,
            double gamma_R,
            double* pressure_hint,
            double& face_speed, double& face_velocity
        ) const
        {
            if (direction < 0 || direction >= hd_core::NDIM)
            {
                throw std::out_of_range(
                    "ExactSolver: invalid direction"
                );
            }

            // Snapshot the global model once: waves and final energy conversion
            // must use the same EoS. No changes to the RiemannSolver interface.
            const auto model = hd_microphysics::eos_type;
            using hd_microphysics::Equation_of_State;
            switch (model)
            {
                case Equation_of_State::IdealGas_Default:
                case Equation_of_State::IdealGas_Mixed:
                case Equation_of_State::Ionized:
                case Equation_of_State::Ionized_Radiative:
                    break;
                default:
                    throw std::invalid_argument("ExactSolver: unsupported equation of state");
            }
            WL.validate(gamma_L);
            WR.validate(gamma_R);
            const bool vacuum_L = (WL[hd_core::RHO] == 0.0);
            const bool vacuum_R = (WR[hd_core::RHO] == 0.0);
            if (vacuum_L && vacuum_R)
            {
                face_speed = 0.0;
                face_velocity = 0.0;
                if (pressure_hint) *pressure_hint = 0.0;
                return hd_core::Conserved{};
            }
            std::optional<PressureWave> left, right;
            if (!vacuum_L) left.emplace(WL, gamma_L, model, tolerance_);
            if (!vacuum_R) right.emplace(WR, gamma_R, model, tolerance_);
            double solved_pressure = 0.0;
            double min_wave_speed = 0.0;
            double max_wave_speed = 0.0;
            InterfaceState interfaceState;
            if(vacuum_R)
            {
                //left gas to right vacuum
                interfaceState = gas_to_vacuum(WL,*left,Side::Left,direction);
                const double u = WL[hd_core::VEL(direction)];
                const double c = left->sound_speed;
                min_wave_speed = u - c;
                max_wave_speed = u - left->vacuum_limit();
            }
            else if(vacuum_L)
            {
                //right gas to left vacuum
                interfaceState = gas_to_vacuum(WR,*right,Side::Right,direction);
                const double u = WR[hd_core::VEL(direction)];
                const double c = right->sound_speed;
                min_wave_speed = u + right->vacuum_limit();
                max_wave_speed = u + c;
            }
            else
            {

                auto star = solve_star(WL,WR,*left,*right,tolerance_,max_iterations_,direction,method_,
                                       pressure_hint ? *pressure_hint : 0.0);
                if(star.has_value())
                {
                    auto [p_star, vn_star] = *star;
                    solved_pressure = p_star;

                    // 无论界面采样在哪一侧，都计算左右两侧波速。
                    min_wave_speed = outer_wave_speed(WL, *left, p_star, Side::Left, direction);
                    max_wave_speed = outer_wave_speed(WR, *right, p_star, Side::Right, direction);

                    interfaceState = sample_interface(
                        WL,WR,
                        *left, *right,
                        vn_star, p_star,direction
                    );
                }
                else
                {
                    // 中间生成真空：最外侧边界仍是左右稀疏波头。
                    min_wave_speed = outer_wave_speed(WL, *left, 0.0, Side::Left, direction);
                    max_wave_speed = outer_wave_speed(WR, *right, 0.0, Side::Right, direction);
                    interfaceState = sample_vacuum(
                        WL,WR,
                        *left,*right,
                        direction
                    );
                }
            }
            if (!std::isfinite(min_wave_speed) ||
                !std::isfinite(max_wave_speed) ||
                min_wave_speed > max_wave_speed)
                throw std::runtime_error("ExactSolver: invalid wave-speed bounds");

            // 即使界面通量为零，也把本次最大传播速度交给网格循环。
            face_velocity = interfaceState.W[hd_core::VEL(direction)];
            face_speed = std::max(std::abs(min_wave_speed), std::abs(max_wave_speed));
            hd_core::Conserved result{};
            if (interfaceState.W[hd_core::RHO] != 0.0)
            {
                // Convert this sampled state once; flux uses the resulting energies directly.
                hd_core::Cell context;
                context.gamma = interfaceState.gamma;
                const auto U = context.to_Conserved(interfaceState.W, model);
                result = physical_flux(interfaceState.W, U, direction);
            }
            if (pressure_hint) *pressure_hint = solved_pressure;
            return result;
        }
}
