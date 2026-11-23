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
#include "microphysics/eos/nuclear.h"
#include "physcon/unit.h"
#include "flux.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>
namespace hydro
{
    namespace
    {
        using hd_microphysics::Equation_of_State;
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
        // One side's ORIGINAL state, unchanged throughout the star-pressure solve.
        // Composition/Saha are frozen. Only the radiation branch needs mu, obtained
        // once from the incoming face state. No composition scan occurs in a wave
        // iteration. Constant specific ionization energy cancels in the Hugoniot.
        struct PressureWave
        {
            double rho0, p0, gamma, sound_speed;
            double shock_a, shock_b;
            double rarefaction_exponent, rarefaction_coefficient;
            Equation_of_State model;
            double gas_fraction = 1.0, radiation_fraction = 0.0;
            double energy_scale, velocity_scale, e0, inner_tolerance;
            mutable std::optional<double> vacuum_value;
            const double rho_n = hd_microphysics::nuclear_density_cgs/physcon::units.unit_density;
            double p_th0;


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
                : rho0(W[hd_core::RHO]), p0(W[hd_core::PRE]), gamma(g),
                  sound_speed(std::sqrt(g * p0 / rho0)),
                  shock_a(2 / ((g + 1) * rho0)),
                  shock_b((g - 1) / (g + 1) * p0),
                  rarefaction_exponent((g - 1) / (2 * g)),
                  rarefaction_coefficient(2 * sound_speed / (g - 1)),
                  model(model),
                  energy_scale(p0/rho0), velocity_scale(std::sqrt(energy_scale)),
                  e0(energy_scale/(g-1)),
                  inner_tolerance(std::clamp(0.01*tolerance, 5e-14, 1e-10))
            {
                if (!std::isfinite(sound_speed) || sound_speed <= 0.0)
                    throw std::invalid_argument("ExactSolver: invalid initial sound speed");
                if (model == Equation_of_State::IdealGas_Default
                    || model == Equation_of_State::IdealGas_Mixed
                    || model == Equation_of_State::Ionized)
                    return;
                if(model == Equation_of_State::Nuclear)
                {
                    double p_c0,u_c0;
                    hd_microphysics::nuclear_cold_state(rho0,gamma,p_c0,u_c0);
                    p_th0 = p0 - p_c0;
                    e0 = u_c0/rho0 + p_th0/(hd_microphysics::nuclear_gamma_thermal - 1.0) *rho0;
                    const double gamma_c = rho0 > rho_n ? hd_microphysics::nuclear_gamma_stiff : gamma;
                    sound_speed = sqrt((gamma_c * p_c0 + hd_microphysics::nuclear_gamma_thermal * p_th0)/rho0);
                    return;
                }

                if (model != Equation_of_State::Ionized_Radiative)
                    throw std::invalid_argument("ExactSolver: unsupported equation of state");
                hd_core::Cell context;
                context.U[hd_core::RHO] = rho0;
                for (int k = 0; k < hd_core::NCOMP; ++k)
                    context.U[hd_core::COMP(k)] = rho0 * W[hd_core::COMP(k)];
                // Implicit Jacobian probes perturb rho and rho*X independently.
                // Check normalization at initialization/accepted steps, not here;
                // the wave still requires finite, positive thermodynamic parameters.
                hd_microphysics::evaluate_cell_density(context);
                const double R = physcon::R_g / context.mean_molecular_weight;
                if (!std::isfinite(R) || R <= 0.0 || !std::isfinite(physcon::a_Rad)
                    || physcon::a_Rad <= 0.0)
                    throw std::invalid_argument("ExactSolver: invalid radiation EoS parameters");

                // Scale by the smaller gas-only / radiation-only temperature.
                // Logarithms avoid forming T^4 or rho*R at an extreme scale.
                const double log_b = std::log(rho0) + std::log(R);
                const double log_d = std::log(physcon::a_Rad/3.0);
                const double log_p = std::log(p0);
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
                double pg, pr;
                sound_speed = velocity_scale * std::sqrt(thermo(1.0, 1.0, pg, pr).second);
                if (!std::isfinite(e0) || !std::isfinite(sound_speed) || sound_speed <= 0.0)
                    throw std::runtime_error("ExactSolver: initial radiation state exceeds numerical range");
            }

            // Dimensionless rho/rho0, P/P0, e/(P0/rho0), c^2/(P0/rho0).
            // e contains thermal+radiation energy only; e_ion(X) cancels exactly.
            // Return {e, c^2}; gas/radiation pressures are ordinary outputs.
            std::pair<double, double> thermo(double r, double p, double& pg, double& pr) const
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
                pg = p*b*t;
                pr = p*d*t*t*t*t;
                const double A = pg + 4.0*pr;
                const double B = pg/(gamma-1.0) + 12.0*pr;
                const double e = (pg/(gamma-1.0) + 3.0*pr)/r;
                const double cs2 = (pg + A*(A/B))/r;
                if (!std::isfinite(e) || !std::isfinite(cs2) || cs2 <= 0.0)
                    throw std::runtime_error("ExactSolver: nonfinite radiation wave state");
                return {e, cs2};
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

            // Radiation rarefaction only: q=log(P/P0), y={log(rho/rho0), f/velocity_scale}.
            std::pair<double, double> integrate(double from, double to,
                                                std::pair<double, double> y) const
            {
                if (to > from || to < -640.0)
                    throw std::runtime_error("ExactSolver: rarefaction exceeds numerical range");
                auto rhs = [&](double q, const std::pair<double, double>& state) {
                    const double r = std::exp(state.first), p = std::exp(q);
                    double pg, pr;
                    const double cs2 = thermo(r, p, pg, pr).second;
                    return std::pair<double, double>{p/r/cs2, p/r/std::sqrt(cs2)};
                };
                auto rk4 = [&](double q, const std::pair<double, double>& state, double h) {
                    const auto a = rhs(q, state);
                    const auto b = rhs(q+0.5*h, {state.first+0.5*h*a.first, state.second+0.5*h*a.second});
                    const auto c = rhs(q+0.5*h, {state.first+0.5*h*b.first, state.second+0.5*h*b.second});
                    const auto d = rhs(q+h, {state.first+h*c.first, state.second+h*c.second});
                    return std::pair<double, double>{
                        state.first+h/6.0*(a.first+2*b.first+2*c.first+d.first),
                        state.second+h/6.0*(a.second+2*b.second+2*c.second+d.second)};
                };
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
                    const double error = std::max(std::abs(fine.first-coarse.first),
                        std::abs(fine.second-coarse.second)/std::max(1.0, std::abs(fine.second))) / 15.0;
                    if (!std::isfinite(error))
                        throw std::runtime_error("ExactSolver: nonfinite rarefaction integration error");
                    if (error <= inner_tolerance)
                    {
                        y.first = fine.first+(fine.first-coarse.first)/15.0;
                        y.second = fine.second+(fine.second-coarse.second)/15.0;
                        q = h == to-q ? to : q+h;
                    }
                    const double factor = error > 0.0
                        ? std::clamp(0.9*std::pow(inner_tolerance/error, 0.2), 0.2, 2.0) : 2.0;
                    h = -std::min(4.0, std::abs(h)*factor);
                }
                return y;
            }

            double vacuum_limit() const
            {
                if (model == Equation_of_State::IdealGas_Default
                    || model == Equation_of_State::IdealGas_Mixed
                    || model == Equation_of_State::Ionized)
                    return -rarefaction_coefficient;
                if (model == Equation_of_State::Ionized_Radiative)
                {
                    if (vacuum_value) return *vacuum_value;
                    std::pair<double, double> y{0.0, 0.0};
                    double q = 0.0;
                    for (int n = 0; n < 160; ++n)
                    {
                        y = integrate(q, q-4.0, y);
                        q -= 4.0;
                        double pg, pr;
                        const double c = std::sqrt(thermo(std::exp(y.first), std::exp(q), pg, pr).second);
                        const double gmin = minimum_gamma1(pg/(pg+pr));
                        // Bound the remaining velocity integral; never query P=0.
                        if (gmin > 1.0 && 2.0*c/(gmin-1.0)
                            <= inner_tolerance*std::max(1.0, sound_speed/velocity_scale))
                        {
                            vacuum_value = velocity_scale*y.second;
                            return *vacuum_value;
                        }
                    }
                    throw std::runtime_error("ExactSolver: vacuum limit did not converge within numerical range");
                }
                throw std::invalid_argument("ExactSolver: unsupported equation of state");
            }

            double radiation_fundamental_derivative() const
            {
                const double beta = gas_fraction/(gas_fraction+radiation_fraction);
                const double A = 4.0-3.0*beta;
                const double B = beta/(gamma-1.0)+12.0*(1.0-beta);
                const double G1 = beta+A*A/B;
                const double dG1_db = 1.0-6.0*A/B
                    - (A/B)*(A/B)*(1.0/(gamma-1.0)-12.0);
                const double db_dlnrho = -beta*(1.0-beta)*(3.0*A/B-1.0);
                return 1.0+0.5*(G1-1.0+dG1_db*db_dlnrho/G1);
            }

            // State on this side's wave curve at a supplied pressure, in code units.
            // This does not solve for star pressure or cache any sampled states.
            // All interface sampling obtains density and sound speed here.
            std::pair<double, double> wave_state_at_pressure(double p) const
            {
                if (!std::isfinite(p) || p < 0.0)
                    throw std::invalid_argument("ExactSolver: invalid wave-state pressure");
                if (model == Equation_of_State::Nuclear)
                {
                    if(p==p0)
                    {
                        return {rho0,sound_speed};
                    }
                    if(p<p0)
                    {
                        
                    }
                    if(p>p0)
                    {

                    }
                }
                if (model != Equation_of_State::IdealGas_Default
                    && model != Equation_of_State::IdealGas_Mixed
                    && model != Equation_of_State::Ionized
                    && model != Equation_of_State::Ionized_Radiative)
                    throw std::invalid_argument("ExactSolver: unsupported equation of state");
                if (p == 0.0) return {0.0, 0.0};
                if (p == p0) return {rho0, sound_speed};
                
                if (model != Equation_of_State::Ionized_Radiative)
                {
                    const double ratio = p/p0;
                    if (p > p0)
                    {
                        const double b = (gamma-1.0)/(gamma+1.0);
                        const double density = rho0*(ratio+b)/(b*ratio+1.0);
                        return {density, std::sqrt(gamma*p/density)};
                    }
                    return {rho0*std::pow(ratio, 1.0/gamma),
                            sound_speed*std::pow(ratio, rarefaction_exponent)};
                }

                const double z = p/p0, dz = (p-p0)/p0;
                if (std::abs(dz) < 1e-7)
                {
                    const double fundamental = radiation_fundamental_derivative();
                    const double x = (p-p0)/(rho0*sound_speed*sound_speed);
                    const double r = 1.0+x-(fundamental-1.0)*x*x;
                    double pg, pr;
                    return {rho0*r, velocity_scale*std::sqrt(thermo(r, z, pg, pr).second)};
                }
                if (p < p0)
                {
                    const auto y = integrate(0.0, std::log(p)-std::log(p0), {0.0, 0.0});
                    const double r = std::exp(y.first);
                    double pg, pr;
                    return {rho0*r, velocity_scale*std::sqrt(thermo(r, z, pg, pr).second)};
                }

                // H(r,z)=e(r,z)-e0-(z+1)/2*(1-1/r)=0, r>=1.
                const double initial_e = e0/energy_scale;
                auto H = [&](double r) {
                    double pg, pr;
                    return thermo(r, z, pg, pr).first - initial_e - 0.5*(z+1.0)*((r-1.0)/r);
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
                const double r = lo+0.5*(hi-lo);
                double pg, pr;
                return {rho0*r, velocity_scale*std::sqrt(thermo(r, z, pg, pr).second)};
            }

            // Velocity jump and pressure derivative only; no sampled-state outputs.
            std::pair<double, double> evaluate_radiation(double p) const
            {
                if (p == 0.0) return {vacuum_limit(), 0.0};
                if (p == p0) return {0.0, 1.0/(rho0*sound_speed)};
                const double z = p/p0, dz = (p-p0)/p0;
                // Weak-wave limit avoids cancellation in the Hugoniot relation.
                if (std::abs(dz) < 1e-7)
                {
                    const double fundamental = radiation_fundamental_derivative();
                    const double x = (p-p0)/(rho0*sound_speed*sound_speed);
                    return {(p-p0)/(rho0*sound_speed)*(1.0-0.5*fundamental*x),
                            (1.0-fundamental*x)/(rho0*sound_speed)};
                }
                const auto [density, c] = wave_state_at_pressure(p);
                if (p < p0)
                {
                    const auto y = integrate(0.0, std::log(p)-std::log(p0), {0.0, 0.0});
                    return {velocity_scale*y.second, 1.0/(density*c)};
                }
                const double r = density/rho0, D = (r-1.0)/r;
                double pg, pr;
                thermo(r, z, pg, pr);
                const double A = pg+4.0*pr, B = pg/(gamma-1.0)+12.0*pr;
                const double e_p = (B/A)/r;
                const double e_rho = -(3.0*pr/r+(B/A)*(pg/r))/r;
                const double f = std::sqrt(dz*D);
                const double dr_dp = (0.5*D-e_p)/(e_rho-0.5*(z+1.0)/r/r);
                const double df_dp = (D+dz*dr_dp/r/r)/(2.0*f);
                if (!std::isfinite(df_dp) || df_dp <= 0.0)
                    throw std::runtime_error("ExactSolver: invalid shock pressure derivative");
                return {velocity_scale*f, velocity_scale/p0*df_dp};
            }

            // Return {wave value, pressure derivative}, sharing sqrt/pow work.
            std::pair<double, double> evaluate(double trial_pressure) const
            {
                if (!std::isfinite(trial_pressure) || trial_pressure < 0.0)
                    throw std::invalid_argument("ExactSolver: invalid trial pressure");
                if (model == Equation_of_State::IdealGas_Default
                    || model == Equation_of_State::IdealGas_Mixed
                    || model == Equation_of_State::Ionized)
                {
                    if (trial_pressure > p0)
                    {
                        const double root = std::sqrt(shock_a / (trial_pressure + shock_b));
                        return {
                            (trial_pressure - p0) * root,
                            root * (1 - (trial_pressure - p0) / (2 * (trial_pressure + shock_b)))
                        };
                    }
                    const double power = std::pow(trial_pressure / p0, rarefaction_exponent);
                    return {
                        rarefaction_coefficient * (power - 1),
                        trial_pressure > 0
                            ? rarefaction_coefficient * rarefaction_exponent * power / trial_pressure
                            : 0
                    };
                }
                if (model == Equation_of_State::Ionized_Radiative)
                {
                    return evaluate_radiation(trial_pressure);
                }
                throw std::invalid_argument("ExactSolver: unsupported equation of state");
            }
        };

        // 非真空气体一侧的最外层波：激波速度或稀疏波头速度。
        /**
         * @brief Compute the outer shock or rarefaction-head speed on one gas side.
         * @param[in] W Non-vacuum primitive state on this side.
         * @param[in] wave Original state and EoS on this side.
         * @param[in] p_star Nonnegative star-region pressure.
         * @param[in] side Left or right side of the Riemann problem.
         * @param[in] direction Normal axis in [0, hd_core::NDIM).
         * @return Signed speed of the outermost wave on this side.
         * @pre W has positive density and pressure and all inputs are physically valid.
         */
        double outer_wave_speed(const hd_core::Primitives& W, const PressureWave& wave,
                                double p_star, Side side, int direction)
        {
            if (wave.model == Equation_of_State::IdealGas_Default
                || wave.model == Equation_of_State::IdealGas_Mixed
                || wave.model == Equation_of_State::Ionized)
            {
                const double gamma = wave.gamma;
                const double p = W[hd_core::PRE];
                const double c = wave.wave_state_at_pressure(p).second;
                const double factor = p_star > p
                    ? std::sqrt((gamma + 1.0) / (2.0 * gamma) * (p_star / p)
                                + (gamma - 1.0) / (2.0 * gamma))
                    : 1.0;
                return W[hd_core::VEL(direction)]
                    + (side == Side::Left ? -1.0 : 1.0) * c * factor;
            }
            if (wave.model == Equation_of_State::Ionized_Radiative)
            {
                double speed = wave.wave_state_at_pressure(wave.p0).second;
                if (p_star > wave.p0)
                    speed = ((p_star-wave.p0)/wave.rho0)/wave.evaluate(p_star).first;
                return W[hd_core::VEL(direction)] + (side == Side::Left ? -speed : speed);
            }
            throw std::invalid_argument("ExactSolver: unsupported equation of state");
        }
        /**
         * @brief Evaluate one side of the exact Riemann pressure equation.
         * @param[in] p_star Trial star pressure, nonnegative.
         * @param[in] wave Original state and EoS on this side.
         * @return Velocity jump contribution from the shock or rarefaction branch.
         * @pre Inputs describe a finite non-vacuum state supported by the EoS.
         */
        double wave_condition(double p_star, const PressureWave& wave);
        /**
         * @brief Solve for star pressure and contact velocity by safeguarded Newton iteration (or reference bisection).
         * @param[in] W_L Original left non-vacuum primitive state.
         * @param[in] W_R Original right non-vacuum primitive state.
         * @param[in] left Original left state and EoS.
         * @param[in] right Original right state and EoS.
         * @param[in] tolerance Finite positive relative pressure-bracket tolerance.
         * @param[in] max_iterations Positive maximum number of pressure iterations.
         * @param[in] direction Normal axis in [0, hd_core::NDIM).
         * @param[in] method Safeguarded Newton or reference bisection.
         * @return Optional pair {star pressure, normal contact velocity}; nullopt when the two gas states generate vacuum.
         * @pre Both input states are valid gas states with positive density and pressure.
         * @throws std::runtime_error If pressure bracketing overflows or bisection does not converge.
         */
        std::optional<std::pair<double, double>> solve_star(
            const hd_core::Primitives& W_L, const hd_core::Primitives& W_R,
            const PressureWave& left, const PressureWave& right, double tolerance,
            int max_iterations, int direction, ExactPressureMethod method, double pressure_hint);
        /**
         * @brief Sample a rarefaction wave at the stationary interface.
         * @param[in] K Original non-vacuum primitive state on the chosen side.
         * @param[in] wave Original state and EoS on this side.
         * @param[in] v_nstar Normal velocity in the star region or at the vacuum edge.
         * @param[in] pre_star Star pressure; zero is allowed for expansion into vacuum.
         * @param[in] side Side of the Riemann problem being sampled.
         * @param[in] direction Normal axis in [0, hd_core::NDIM).
         * @return Primitive state at similarity coordinate x/t = 0; tangential velocities are inherited from K.
         * @pre Inputs satisfy the selected wave branch; K has positive density and pressure.
         */
        hd_core::Primitives sample_rarefraction(const hd_core::Primitives& K, const PressureWave& wave, double v_nstar, double pre_star,Side side, int direction);
        /**
         * @brief Sample a shock wave at the stationary interface.
         * @param[in] K Original non-vacuum primitive state on the chosen side.
         * @param[in] wave Original state and EoS on this side.
         * @param[in] v_nstar Normal velocity in the star region or at the vacuum edge.
         * @param[in] pre_star Star pressure greater than the original gas pressure.
         * @param[in] side Side of the Riemann problem being sampled.
         * @param[in] direction Normal axis in [0, hd_core::NDIM).
         * @return Primitive state at similarity coordinate x/t = 0; tangential velocities are inherited from K.
         * @pre Inputs satisfy the selected wave branch; K has positive density and pressure.
         */
        hd_core::Primitives sample_shock(const hd_core::Primitives& K, const PressureWave& wave, double v_nstar, double pre_star, Side side, int direction);
        /**
         * @brief Choose the contact side and sample the exact gas interface state.
         * @param[in] L Original left non-vacuum primitive state.
         * @param[in] R Original right non-vacuum primitive state.
         * @param[in] left Original left state and EoS.
         * @param[in] right Original right state and EoS.
         * @param[in] v_nstar Normal contact velocity.
         * @param[in] pre_star Positive converged star pressure.
         * @param[in] direction Normal axis in [0, hd_core::NDIM).
         * @return Sampled primitive state and the adiabatic index of its originating side.
         * @pre Both states and the solved star quantities are valid.
         * @note A strictly positive contact velocity selects the left side; zero selects the right.
         */
        InterfaceState sample_interface(const hd_core::Primitives& L, const hd_core::Primitives& R, const PressureWave& left, const PressureWave& right, double v_nstar, double pre_star, int direction);
        double wave_condition(double p_star, const PressureWave& wave)
        {
            if (wave.model == Equation_of_State::IdealGas_Default
                || wave.model == Equation_of_State::IdealGas_Mixed
                || wave.model == Equation_of_State::Ionized)
            {
                const double rho_K = wave.rho0, p_K = wave.p0, gamma_K = wave.gamma;
                if(p_star>p_K)
                {
                    //Shock
                    double A_K = 2/((gamma_K+1)*rho_K);
                    double B_K = (gamma_K-1)/(gamma_K+1)*p_K;
                    return (p_star-p_K)*std::sqrt(A_K/(p_star+B_K));
                }
                else
                {
                    //Rarefraction
                    double c_K = wave.wave_state_at_pressure(p_K).second;
                    return 2*c_K/(gamma_K-1)*(std::pow(p_star/p_K,(gamma_K-1)/(2*gamma_K))-1);
                }
            }
            if (wave.model == Equation_of_State::Ionized_Radiative)
            {
                return wave.evaluate(p_star).first;
            }
            throw std::invalid_argument("ExactSolver: unsupported equation of state");
        }
        hd_core::Primitives sample_rarefraction(const hd_core::Primitives& K, const PressureWave& wave,
                                                double v_nstar, double pre_star, Side side, int direction)
        {
            if (wave.model == Equation_of_State::IdealGas_Default
                || wave.model == Equation_of_State::IdealGas_Mixed
                || wave.model == Equation_of_State::Ionized)
            {
                const double gamma_K = wave.gamma;
                const double g = gamma_K;
                auto v_nK = K[hd_core::VEL(direction)];
                auto pre_K = K[hd_core::PRE];
                auto c_sK = wave.wave_state_at_pressure(pre_K).second;
                const auto [rho_star, c_sK_star] = wave.wave_state_at_pressure(pre_star);
                auto v_rare_head = v_nK + (side==Side::Left? -1.0 : 1.0) * c_sK;
                auto v_rare_tail = v_nstar + (side==Side::Left? -1.0 : 1.0) * c_sK_star;
                if((side==Side::Left && v_rare_head >=0) || (side==Side::Right && v_rare_head<=0))
                {
                    return K;
                }
                else
                {
                    if( (side==Side::Left && v_rare_tail <=0 )|| (side==Side::Right && v_rare_tail >=0) )
                    {
                        hd_core::Primitives face = K;
                        face[hd_core::RHO] = rho_star;
                        face[hd_core::VEL(direction)] = v_nstar;
                        face[hd_core::PRE] = pre_star;

                        return face;
                    }
                    else
                    {
                        auto v_nf = 1/(g+1)*((g-1)*v_nK +(side==Side::Left? 1.0:-1.0)* 2 *c_sK);
                        auto pre_f = pre_K * std::pow(
                            ((side == Side::Left) ? v_nf : -v_nf)/c_sK, 2*g/(g-1));
                        hd_core::Primitives face = K;
                        face[hd_core::RHO] = wave.wave_state_at_pressure(pre_f).first;
                        face[hd_core::VEL(direction)] = v_nf;
                        face[hd_core::PRE] = pre_f;
                        return face;
                    }
                }
            }
            if (wave.model == Equation_of_State::Ionized_Radiative)
            {
                const double sign = side == Side::Left ? -1.0 : 1.0;
                const double u0 = K[hd_core::VEL(direction)];
                const double c0 = wave.wave_state_at_pressure(wave.p0).second;
                if (sign*(u0+sign*c0) <= 0.0) return K;
                hd_core::Primitives face = K;
                const auto [rho_star, c_star] = wave.wave_state_at_pressure(pre_star);
                if (sign*(v_nstar + sign*c_star) >= 0.0)
                {
                    face[hd_core::RHO] = rho_star;
                    face[hd_core::VEL(direction)] = v_nstar;
                    face[hd_core::PRE] = pre_star;
                    return face;
                }
                // The face lies inside the fan. Solve u(P)+sign*c(P)=x/t=0.
                // Use log pressure so near-vacuum fans do not need linear-P bisection.
                auto characteristic = [&](double q) {
                    const double p = wave.p0*std::exp(q);
                    const double c = wave.wave_state_at_pressure(p).second;
                    const auto value = wave.evaluate(p);
                    return sign*u0 + value.first + c;
                };
                double lo = pre_star > 0.0 ? std::log(pre_star)-std::log(wave.p0) : -1.0;
                double hi = 0.0;
                if (pre_star == 0.0)
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
                const double p = wave.p0*std::exp(lo+0.5*(hi-lo));
                face[hd_core::RHO] = wave.wave_state_at_pressure(p).first;
                const auto value = wave.evaluate(p);
                face[hd_core::VEL(direction)] = u0 + sign*value.first;
                face[hd_core::PRE] = p;
                return face;
            }
            throw std::invalid_argument("ExactSolver: unsupported equation of state");
        }
        hd_core::Primitives sample_shock(const hd_core::Primitives& K, const PressureWave& wave,
                                         double v_nstar, double pre_star, Side side, int direction)
        {
            if (wave.model == Equation_of_State::IdealGas_Default
                || wave.model == Equation_of_State::IdealGas_Mixed
                || wave.model == Equation_of_State::Ionized)
            {
                const double gamma_K = wave.gamma;
                const double g = gamma_K;
                auto v_nK = K[hd_core::VEL(direction)];
                auto pre_K = K[hd_core::PRE];
                auto comp_ratio = pre_star / pre_K;
                auto c_sK = wave.wave_state_at_pressure(pre_K).second;
                auto v_shock = v_nK + (side==Side::Left? -1.0 : 1.0) * c_sK * std::sqrt((g+1)/(2*g)*comp_ratio+(g-1)/(2*g));
                if((v_shock>0 && side==Side::Left) || (v_shock<0 && side==Side::Right))
                {
                    return K;
                }
                else
                {
                    hd_core::Primitives face = K;
                    face[hd_core::RHO] = wave.wave_state_at_pressure(pre_star).first;
                    face[hd_core::VEL(direction)] = v_nstar;
                    face[hd_core::PRE] = pre_star;
                    return face;
                }
            }
            if (wave.model == Equation_of_State::Ionized_Radiative)
            {
                const double speed = outer_wave_speed(K, wave, pre_star, side, direction);
                if ((side == Side::Left && speed > 0.0) || (side == Side::Right && speed < 0.0)) return K;
                hd_core::Primitives face = K;
                face[hd_core::RHO] = wave.wave_state_at_pressure(pre_star).first;
                face[hd_core::VEL(direction)] = v_nstar;
                face[hd_core::PRE] = pre_star;
                return face;
            }
            throw std::invalid_argument("ExactSolver: unsupported equation of state");
        }
        InterfaceState sample_interface(const hd_core::Primitives& L, const hd_core::Primitives& R, const PressureWave& left, const PressureWave& right, double v_nstar, double pre_star, int direction)
        {
            const double gamma_L = left.gamma, gamma_R = right.gamma;
            auto pre_L = L[hd_core::PRE], pre_R = R[hd_core::PRE];
            if(v_nstar>0)
            {
                const double comp_ratio = pre_star / pre_L;
                if(comp_ratio > 1)
                {
                    return {sample_shock(L,left,v_nstar,pre_star,Side::Left,direction),gamma_L};
                }
                else
                {
                    return {sample_rarefraction(L,left,v_nstar,pre_star,Side::Left,direction),gamma_L};
                }
            }
            else
            {
                const double comp_ratio = pre_star / pre_R;
                if(comp_ratio>1)
                {
                    return {sample_shock(R,right,v_nstar,pre_star,Side::Right,direction),gamma_R};
                }
                else
                {
                    return {sample_rarefraction(R,right,v_nstar,pre_star,Side::Right,direction),gamma_R};
                }
            }
        }
        //Original Version, not simplified...
        // hd_core::Primitive sample_interface(const hd_core::Primitive& L, const hd_core::Primitive& R, double gamma_L,  double gamma_R, double v_nstar, double pre_star, Side side)
        // {
        //     auto rho_L = L.density, rho_R = R.density;
        //     auto v_nL = L.velocity_normal, v_nR = R.velocity_normal;
        //     auto pre_L = L.pressure, pre_R = R.pressure;
        //     auto c_sL = std::sqrt(gamma_L * pre_L/rho_L);
        //     auto c_sR = std::sqrt(gamma_R * pre_R/rho_R);
        //     if(v_nstar>0)
        //     {

        //         ///The grid interface is left to the contact interface
        //         const double g = gamma_L;
        //         const double b = (g - 1) / (g + 1);
        //         const double comp_ratio = pre_star / pre_L;
        //         if(comp_ratio > 1)
        //         {
        //             //left shock(left to the contact)
        //             auto v_shock = v_nL - c_sL*std::sqrt((g+1)/(2*g)*comp_ratio+(g-1)/(2*g));
        //             if(v_shock>0)
        //             {
        //             //left shock is located to the right of grid interface(positive)

        //             //  Original Left  grid interface(0) left shock   left star      contact   right star      rarefraction tail    head   original right
        //             //----------------------X----------------|***********************|||*************************|~~~~~~~~~~~~~~~~~~~|-------------------------
        //                 return L;
        //             }
        //             else
        //             {
        //             //left shock is located to the left of grid interface(negative)----interface locates inside the star region(shock-propagation-side)

        //             //  Original Left  left shock   grid interface(0)  left star      contact   right star      rarefraction tail    head   original right
        //             //----------------------|***************X***********************|||*************************|~~~~~~~~~~~~~~~~~~~|-------------------------
        //                 auto rho_star = rho_L * (comp_ratio + b) / ( b * comp_ratio + 1);
        //                 return {rho_star, v_nstar, pre_star};
        //             }
        //         }
        //         else
        //         {
        //             //Remember: the rarefraction starts from the **CONTACT INTERFACE**, so the head of rarefraction
        //             //is the side that's FURTHER away from contact interface, and the tail is the closer side
        //             //left rarefraction(left to the contact)
        //             auto c_sL_star = c_sL * std::pow(pre_star/pre_L,(g-1)/(2*g));
        //             auto v_rare_head = v_nL - c_sL;
        //             auto v_rare_tail = v_nstar - c_sL_star;
        //             if(v_rare_head>=0)
        //             {
        //             //left rarefraction is located to the right of grid interface(positive)

        //             //  Original Left  grid interface(0) rarefraction head     tail  left star    contact    right star       right shock   original right
        //             //----------------------X----------------|~~~~~~~~~~~~~~~~~~~|*****************|||***************************|-----------------------------
        //             //interface at original left
        //                 return L;
        //             }
        //             else
        //             {
        //                 if(v_rare_tail<=0)
        //                 {
        //                 //lThe whole rarefraction is right to the grid interface, so the interface locates inside the star region(rarefraction-propagation-side)

        //                 //  Original Left     rf head              rf tail left star grid interface     contact    right star       right shock   original right
        //                 //----------------------|~~~~~~~~~~~~~~~~~~~|*********************X**************|||****************************|-----------------------------
        //                     auto rho_star = rho_L * std::pow(pre_star / pre_L, 1.0 / g);
        //                     return {rho_star,v_nstar,pre_star};
        //                 }
        //                 else
        //                 {
        //                     auto v_nf = 1/(g+1)*((g-1)*v_nL + 2*c_sL);
        //                     auto c_sf = v_nf;
        //                     auto rho_f = rho_L * std::pow((c_sf/c_sL),2 / (g - 1));
        //                     auto pre_f = pre_L * std::pow((c_sf/c_sL),2 * g / (g-1));
        //                     return {rho_f,v_nf,pre_f};
        //                 }
        //             }
        //         }
        //     }
        //     else
        //     {
        //         const double g = gamma_R;
        //         const double b = (g - 1) / (g + 1);
        //         const double comp_ratio = pre_star / pre_R;
        //         if(comp_ratio > 1)
        //         {
        //             //right shock(right to the contact)
        //             auto v_shock = v_nR + c_sR*std::sqrt((g+1)/(2*g)*comp_ratio+(g-1)/(2*g));
        //             if(v_shock<0)
        //             {
        //             //right shock is located to the left of grid interface(negative)

        //             //-------------|~~~~~~~~~~~~~~|*********||||*******|-------X-------
        //                 return R;
        //             }
        //             else
        //             {
        //             //right shock is located to the right of grid interface, while contact in the left
        //             //Right shock-side star region
        //             //--------------|~~~~~~~~~~~~|*********|||*******X********|-------------------
        //                 auto rho_star = rho_R * (comp_ratio + b) / ( b * comp_ratio + 1);
        //                 return {rho_star, v_nstar, pre_star};
        //             }
        //         }
        //         else
        //         {
        //             //Remember: the rarefraction starts from the **CONTACT INTERFACE**, so the head of rarefraction
        //             //is the side that's FURTHER away from contact interface, and the tail is the closer side
        //             //right rarefraction(right to the contact)
        //             auto c_sR_star = c_sR * std::pow(pre_star/pre_R,(g-1)/(2*g));
        //             auto v_rare_head = v_nR + c_sR;
        //             auto v_rare_tail = v_nstar + c_sR_star;
        //             if(v_rare_head<=0)
        //             {
        //             //right rarefraction is located to the left of grid interface(positive)
        //             //-------|***********|||-----------|~~~~~~~~~~~~~~|---------X------
        //                 return R;
        //             }
        //             else
        //             {
        //                 if(v_rare_tail>=0)
        //                 {
        //                 //The whole rarefraction is left to the grid interface, so the interface locates inside the star region(rarefraction-propagation-side)

        //                 //-----|*******|||*****X*****|~~~~~~~|--------
        //                     auto rho_star = rho_R * std::pow(pre_star / pre_R, 1.0 / g);
        //                     return {rho_star,v_nstar,pre_star};
        //                 }
        //                 else
        //                 {
        //                     //------|*******|||******|~~~~X~~~~~~|-------
        //                     auto v_nf = 1/(g+1)*((g-1)*v_nR - 2*c_sR);
        //                     auto c_sf = -v_nf;
        //                     auto rho_f = rho_R * std::pow((c_sf/c_sR),2 / (g - 1));
        //                     auto pre_f = pre_R * std::pow((c_sf/c_sR),2 * g / (g-1));
        //                     return {rho_f,v_nf,pre_f};
        //                 }
        //             }
        //         }
        //     }

        // }

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
            if (left.model == Equation_of_State::Ionized_Radiative
                && pre_L == pre_R && v_nL == v_nR)
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

        // Keep the reference arithmetic and evaluation order for byte reproducibility.
        // This separate loop needs neither derivative evaluation nor Newton safeguards.
        std::optional<std::pair<double, double>> solve_star_bisection(
            const hd_core::Primitives& W_L, const hd_core::Primitives& W_R,
            const PressureWave& left, const PressureWave& right, double tolerance,
            int max_iterations, int direction)
        {
            RootSolveRecord record;
            const double pre_L = W_L[hd_core::PRE], pre_R = W_R[hd_core::PRE];
            const double v_nL = W_L[hd_core::VEL(direction)], v_nR = W_R[hd_core::VEL(direction)];
            if (left.model == Equation_of_State::Ionized_Radiative
                && pre_L == pre_R && v_nL == v_nR)
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

        /**
         * @brief Sample a vacuum-producing Riemann solution at the stationary interface.
         * @param[in] W_L Original left gas state with positive density and pressure.
         * @param[in] W_R Original right gas state with positive density and pressure.
         * @param[in] left Original left state and EoS.
         * @param[in] right Original right state and EoS.
         * @param[in] direction Normal axis in [0, hd_core::NDIM).
         * @return Vacuum or a sampled rarefaction state, together with its adiabatic index.
         * @pre The two valid gas states satisfy the vacuum-generation criterion.
         */
        InterfaceState sample_vacuum(const hd_core::Primitives& W_L, const hd_core::Primitives& W_R,
                                     const PressureWave& left, const PressureWave& right, int direction)
        {
            if (left.model == Equation_of_State::IdealGas_Default
                || left.model == Equation_of_State::IdealGas_Mixed
                || left.model == Equation_of_State::Ionized)
            {
                const double gamma_L = left.gamma, gamma_R = right.gamma;
                //deal with vacuum generated in Riemann Problem
                //-----------|~~~~~~~~|   vacuum    |~~~~~~~~|-------------
                const double v_nL = W_L[hd_core::VEL(direction)];
                const double v_nR = W_R[hd_core::VEL(direction)];
                const double c_sL = left.wave_state_at_pressure(left.p0).second;
                const double c_sR = right.wave_state_at_pressure(right.p0).second;
                const double S_VL = v_nL + 2.0 * c_sL / (gamma_L - 1.0);//(Left-Vacuum Interace velocity)
                const double S_VR = v_nR - 2.0 * c_sR / (gamma_R - 1.0);
                if(S_VL<=0 && S_VR>=0)
                {
                    return {hd_core::Primitives{},gamma_L}; // In the vacuum
                }
                if(S_VL>0.0)
                {
                    return {sample_rarefraction(W_L,left,S_VL,0.0,Side::Left,direction),gamma_L};
                }
                return {sample_rarefraction(W_R,right,S_VR,0.0,Side::Right,direction),gamma_R};
                //Though there is no star region. In this case, the interface between the tail of rarefraction and vacuum is equivalent to star region, because they are both the boundary of rarefraction(tail).
                //So S_V is v_nstar because it's the speed of boundary, and p_star =0.0, because it can make c_sK_star = 0, the v_rare_tail = S_V....everything fine!
            }
            if (left.model == Equation_of_State::Ionized_Radiative)
            {
                const double S_VL = W_L[hd_core::VEL(direction)]-left.vacuum_limit();
                const double S_VR = W_R[hd_core::VEL(direction)]+right.vacuum_limit();
                if (S_VL <= 0.0 && S_VR >= 0.0) return {hd_core::Primitives{}, left.gamma};
                if (S_VL > 0.0)
                    return {sample_rarefraction(W_L,left,S_VL,0.0,Side::Left,direction),left.gamma};
                return {sample_rarefraction(W_R,right,S_VR,0.0,Side::Right,direction),right.gamma};
            }
            throw std::invalid_argument("ExactSolver: unsupported equation of state");
        }
        /**
         * @brief Sample expansion of a gas state into an initially empty neighbor.
         * @param[in] gas Non-vacuum primitive gas state.
         * @param[in] wave Original gas state and EoS.
         * @param[in] gas_side Side initially occupied by gas.
         * @param[in] direction Normal axis in [0, hd_core::NDIM).
         * @return Vacuum or a rarefaction state at x/t = 0, with the gas adiabatic index.
         * @pre gas has finite positive density and pressure.
         */
        InterfaceState gas_to_vacuum(const hd_core::Primitives& gas, const PressureWave& wave,
                                     Side gas_side, int direction)
        {
            if (wave.model == Equation_of_State::IdealGas_Default
                || wave.model == Equation_of_State::IdealGas_Mixed
                || wave.model == Equation_of_State::Ionized)
            {
                const double gamma = wave.gamma;
                //To deal with real vacuum cell
                //Basic Idea(take left gas, right vacuum as example)
                //Gas expand to vacuum, producing rarefraction
                //--------------------
                const double v_n = gas[hd_core::VEL(direction)];
                const double c_s = wave.wave_state_at_pressure(wave.p0).second;
                const double S_V = v_n + (gas_side == Side::Left? 1.0 : -1.0) * 2.0 * c_s / (gamma - 1.0);
                if((gas_side == Side::Left && S_V<= 0.0) || (gas_side==Side::Right && S_V>=0.0))
                {
                    return {hd_core::Primitives{},gamma};
                }
                return {sample_rarefraction(gas,wave,S_V,0.0,gas_side,direction),gamma};
                //Though there is no star region. In this case, the interface between the tail of rarefraction and vacuum is equivalent to star region, because they are both the boundary of rarefraction(tail).
                //So S_V is v_nstar because it's the speed of boundary, and p_star =0.0, because it can make c_sK_star = 0, the v_rare_tail = S_V....everything fine!
            }
            if (wave.model == Equation_of_State::Ionized_Radiative)
            {
                const double S_V = gas[hd_core::VEL(direction)]
                    + (gas_side == Side::Left ? -1.0 : 1.0)*wave.vacuum_limit();
                if ((gas_side == Side::Left && S_V <= 0.0) || (gas_side == Side::Right && S_V >= 0.0))
                    return {hd_core::Primitives{},wave.gamma};
                return {sample_rarefraction(gas,wave,S_V,0.0,gas_side,direction),wave.gamma};
            }
            throw std::invalid_argument("ExactSolver: unsupported equation of state");
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
                const double c = left->wave_state_at_pressure(left->p0).second;
                min_wave_speed = u - c;
                max_wave_speed = u - left->vacuum_limit();
            }
            else if(vacuum_L)
            {
                //right gas to left vacuum
                interfaceState = gas_to_vacuum(WR,*right,Side::Right,direction);
                const double u = WR[hd_core::VEL(direction)];
                const double c = right->wave_state_at_pressure(right->p0).second;
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
