/**
 * @file exact.cpp
 * @brief Hydrodynamic spatial operators.
 */
#include "exact.h"
#include "utils/performance.h"
#include "hydro/wave_timestep.h"
#include "core.h"
#include "flux.h"
#include <algorithm>
#include <cmath>
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
        // 非真空气体一侧的最外层波：激波速度或稀疏波头速度。
        /**
         * @brief Compute the outer shock or rarefaction-head speed on one gas side.
         * @param[in] W Non-vacuum primitive state on this side.
         * @param[in] gamma Adiabatic index greater than one.
         * @param[in] p_star Nonnegative star-region pressure.
         * @param[in] side Left or right side of the Riemann problem.
         * @param[in] direction Normal axis in [0, hd_core::NDIM).
         * @return Signed speed of the outermost wave on this side.
         * @pre W has positive density and pressure and all inputs are physically valid.
         */
        double outer_wave_speed(const hd_core::Primitives& W, double gamma,
                                double p_star, Side side, int direction)
        {
            const double p = W[hd_core::PRE];
            const double c = std::sqrt(gamma * p / W[hd_core::RHO]);
            const double factor = p_star > p
                ? std::sqrt((gamma + 1.0) / (2.0 * gamma) * (p_star / p)
                            + (gamma - 1.0) / (2.0 * gamma))
                : 1.0;
            return W[hd_core::VEL(direction)]
                + (side == Side::Left ? -1.0 : 1.0) * c * factor;
        }
        /**
         * @brief Evaluate one side of the exact Riemann pressure equation.
         * @param[in] p_star Trial star pressure, nonnegative.
         * @param[in] rho_K Positive density of the original gas state.
         * @param[in] p_K Positive pressure of the original gas state.
         * @param[in] gamma_K Adiabatic index greater than one.
         * @return Velocity jump contribution from the shock or rarefaction branch.
         * @pre Inputs describe a finite non-vacuum ideal-gas state.
         */
        double wave_condition(double p_star, double rho_K, double p_K, double gamma_K);
        /**
         * @brief Solve for star pressure and contact velocity by safeguarded Newton iteration (or reference bisection).
         * @param[in] W_L Original left non-vacuum primitive state.
         * @param[in] W_R Original right non-vacuum primitive state.
         * @param[in] gamma_L Left adiabatic index greater than one.
         * @param[in] gamma_R Right adiabatic index greater than one.
         * @param[in] tolerance Finite positive relative pressure-bracket tolerance.
         * @param[in] max_iterations Positive maximum number of pressure iterations.
         * @param[in] direction Normal axis in [0, hd_core::NDIM).
         * @param[in] method Safeguarded Newton or reference bisection.
         * @return Optional pair {star pressure, normal contact velocity}; nullopt when the two gas states generate vacuum.
         * @pre Both input states are valid gas states with positive density and pressure.
         * @throws std::runtime_error If pressure bracketing overflows or bisection does not converge.
         */
        std::optional<std::pair<double, double>> solve_star(const hd_core::Primitives& W_L, const hd_core::Primitives& W_R,double gamma_L, double gamma_R,double tolerance,int max_iterations,int direction, ExactPressureMethod method);
        /**
         * @brief Sample a rarefaction wave at the stationary interface.
         * @param[in] K Original non-vacuum primitive state on the chosen side.
         * @param[in] gamma_K Adiabatic index of K, greater than one.
         * @param[in] v_nstar Normal velocity in the star region or at the vacuum edge.
         * @param[in] pre_star Star pressure; zero is allowed for expansion into vacuum.
         * @param[in] side Side of the Riemann problem being sampled.
         * @param[in] direction Normal axis in [0, hd_core::NDIM).
         * @return Primitive state at similarity coordinate x/t = 0; tangential velocities are inherited from K.
         * @pre Inputs satisfy the selected wave branch; K has positive density and pressure.
         */
        hd_core::Primitives sample_rarefraction(const hd_core::Primitives& K, double gamma_K, double v_nstar, double pre_star,Side side, int direction);
        /**
         * @brief Sample a shock wave at the stationary interface.
         * @param[in] K Original non-vacuum primitive state on the chosen side.
         * @param[in] gamma_K Adiabatic index of K, greater than one.
         * @param[in] v_nstar Normal velocity in the star region or at the vacuum edge.
         * @param[in] pre_star Star pressure greater than the original gas pressure.
         * @param[in] side Side of the Riemann problem being sampled.
         * @param[in] direction Normal axis in [0, hd_core::NDIM).
         * @return Primitive state at similarity coordinate x/t = 0; tangential velocities are inherited from K.
         * @pre Inputs satisfy the selected wave branch; K has positive density and pressure.
         */
        hd_core::Primitives sample_shock(const hd_core::Primitives& K, double gamma_K, double v_nstar, double pre_star, Side side, int direction);
        /**
         * @brief Choose the contact side and sample the exact gas interface state.
         * @param[in] L Original left non-vacuum primitive state.
         * @param[in] R Original right non-vacuum primitive state.
         * @param[in] gamma_L Left adiabatic index greater than one.
         * @param[in] gamma_R Right adiabatic index greater than one.
         * @param[in] v_nstar Normal contact velocity.
         * @param[in] pre_star Positive converged star pressure.
         * @param[in] direction Normal axis in [0, hd_core::NDIM).
         * @return Sampled primitive state and the adiabatic index of its originating side.
         * @pre Both states and the solved star quantities are valid.
         * @note A strictly positive contact velocity selects the left side; zero selects the right.
         */
        InterfaceState sample_interface(const hd_core::Primitives& L, const hd_core::Primitives& R, double gamma_L, double gamma_R, double v_nstar, double pre_star, int direction);
        double wave_condition(double p_star, double rho_K, double p_K, double gamma_K)
        {
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
                double c_K = std::sqrt(gamma_K*p_K/rho_K);
                return 2*c_K/(gamma_K-1)*(std::pow(p_star/p_K,(gamma_K-1)/(2*gamma_K))-1);
            }
        }
        hd_core::Primitives sample_rarefraction(const hd_core::Primitives& K, double gamma_K, double v_nstar, double pre_star, Side side, int direction)
        {
            const double g = gamma_K;
            auto rho_K = K[hd_core::RHO];
            auto v_nK = K[hd_core::VEL(direction)];
            auto pre_K = K[hd_core::PRE];
            auto c_sK = std::sqrt(gamma_K * pre_K/rho_K);
            auto c_sK_star = c_sK * std::pow(pre_star/pre_K,(g-1)/(2*g));
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
                    auto rho_star = rho_K * std::pow(pre_star / pre_K, 1.0 / g);
                    hd_core::Primitives face = K;
                    face[hd_core::RHO] = rho_star;
                    face[hd_core::VEL(direction)] = v_nstar;
                    face[hd_core::PRE] = pre_star;

                    return face;
                }
                else
                {
                    auto v_nf = 1/(g+1)*((g-1)*v_nK +(side==Side::Left? 1.0:-1.0)* 2 *c_sK);
                    auto c_sf = (side == Side::Left) ? v_nf : -v_nf;
                    auto rho_f = rho_K * std::pow((c_sf/c_sK),2 / (g - 1));
                    auto pre_f = pre_K * std::pow((c_sf/c_sK),2 * g / (g-1));
                    hd_core::Primitives face = K;
                    face[hd_core::RHO] = rho_f;
                    face[hd_core::VEL(direction)] = v_nf;
                    face[hd_core::PRE] = pre_f;
                    return face;
                }
            }
        }
        hd_core::Primitives sample_shock(const hd_core::Primitives& K, double gamma_K, double v_nstar, double pre_star, Side side, int direction)
        {
            const double g = gamma_K;
            const double b = (g-1)/(g+1);
            auto rho_K = K[hd_core::RHO];
            auto v_nK = K[hd_core::VEL(direction)];
            auto pre_K = K[hd_core::PRE];
            auto comp_ratio = pre_star / pre_K;
            auto c_sK = std::sqrt(gamma_K * pre_K/rho_K);
            auto v_shock = v_nK + (side==Side::Left? -1.0 : 1.0) * c_sK * std::sqrt((g+1)/(2*g)*comp_ratio+(g-1)/(2*g));
            if((v_shock>0 && side==Side::Left) || (v_shock<0 && side==Side::Right))
            {
                return K;
            }
            else
            {
                hd_core::Primitives face = K;
                auto rho_star = rho_K * (comp_ratio + b) / (b*comp_ratio + 1);
                face[hd_core::RHO] = rho_star;
                face[hd_core::VEL(direction)] = v_nstar;
                face[hd_core::PRE] = pre_star;
                return face;
            }
            
        }
        InterfaceState sample_interface(const hd_core::Primitives& L, const hd_core::Primitives& R, double gamma_L,  double gamma_R, double v_nstar, double pre_star, int direction)
        {
            auto pre_L = L[hd_core::PRE], pre_R = R[hd_core::PRE];
            if(v_nstar>0)
            {
                const double comp_ratio = pre_star / pre_L;
                if(comp_ratio > 1)
                {
                    return {sample_shock(L,gamma_L,v_nstar,pre_star,Side::Left,direction),gamma_L};
                }
                else
                {
                    return {sample_rarefraction(L,gamma_L,v_nstar,pre_star,Side::Left,direction),gamma_L};
                }
            }
            else
            {
                const double comp_ratio = pre_star / pre_R;
                if(comp_ratio>1)
                {
                    return {sample_shock(R,gamma_R,v_nstar,pre_star,Side::Right,direction),gamma_R};
                }
                else
                {
                    return {sample_rarefraction(R,gamma_R,v_nstar,pre_star,Side::Right,direction),gamma_R};
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
        
        std::optional<std::pair<double, double>> solve_star(const hd_core::Primitives& W_L, const hd_core::Primitives& W_R,double gamma_L, double gamma_R,double tolerance,int max_iterations,int direction, ExactPressureMethod method)
        {
            // Aggregate once per solve, including exception paths. No clock or allocation.
            struct Record
            {
                std::uint64_t iterations = 0, evaluations = 0, initial = 0, newton = 0, bisect = 0, expansions = 0;
                bool failed = true, vacuum = false;
                ~Record()
                {
                    if (!performance::enabled()) return;
                    auto& r = performance::exact_roots;
                    ++r.solves; r.iterations += iterations; r.evaluations += evaluations;
                    r.initial_guesses += initial; r.newton_steps += newton; r.bisection_steps += bisect;
                    r.bracket_expansions += expansions; r.failures += failed; r.vacuum += vacuum;
                    r.max_iterations = std::max(r.max_iterations, iterations);
                    const auto bin = iterations == 0 ? 0 : iterations <= 4 ? 1 : iterations <= 8 ? 2 :
                        iterations <= 16 ? 3 : iterations <= 32 ? 4 : iterations <= 64 ? 5 : 6;
                    ++r.histogram[bin];
                }
            } record;
            struct Wave
            {
                double p, c, a, b, exponent, coefficient;
                Wave(double rho, double pressure, double gamma)
                    : p(pressure), c(std::sqrt(gamma * p / rho)),
                      a(2 / ((gamma + 1) * rho)), b((gamma - 1) / (gamma + 1) * p),
                      exponent((gamma - 1) / (2 * gamma)), coefficient(2 * c / (gamma - 1)) {}
                // Value and analytic derivative share sqrt/pow work.
                std::pair<double, double> evaluate(double x) const
                {
                    if (x > p)
                    {
                        const double root = std::sqrt(a / (x + b));
                        return {(x - p) * root, root * (1 - (x - p) / (2 * (x + b)))};
                    }
                    const double power = std::pow(x / p, exponent);
                    return {coefficient * (power - 1), x > 0 ? coefficient * exponent * power / x : 0};
                }
            };
            const double rho_L = W_L[hd_core::RHO], rho_R = W_R[hd_core::RHO];
            const double pre_L = W_L[hd_core::PRE], pre_R = W_R[hd_core::PRE];
            const double v_nL = W_L[hd_core::VEL(direction)], v_nR = W_R[hd_core::VEL(direction)];
            const Wave left(rho_L, pre_L, gamma_L), right(rho_R, pre_R, gamma_R);
            const bool legacy = method == ExactPressureMethod::Bisection;
            auto residual = [&](double p) {
                ++record.evaluations;
                if (legacy) return std::pair<double, double>{
                    wave_condition(p, rho_L, pre_L, gamma_L) +
                    wave_condition(p, rho_R, pre_R, gamma_R) + v_nR - v_nL, 0};
                const auto l = left.evaluate(p), r = right.evaluate(p);
                return std::pair<double, double>{l.first + r.first + v_nR - v_nL, l.second + r.second};
            };
            if (residual(0).first >= 0)
            {
                record.failed = false; record.vacuum = true;
                return std::nullopt;
            }
            double low = 0, high = std::max(pre_L, pre_R);
            while (residual(high).first < 0)
            {
                high *= 2;
                ++record.expansions;
                if (!std::isfinite(high)) throw std::runtime_error("ExactSolver: pressure bracket overflow");
            }
            double p = (low + high) / 2;
            if (!legacy)
            {
                // Acoustic (PVRS) estimate is inexpensive and accurate for smooth flow.
                // It is only a starting guess; all iterates still use the exact equation.
                const double guess = 0.5 * (pre_L + pre_R) - 0.125 * (v_nR - v_nL) *
                    (rho_L + rho_R) * (left.c + right.c);
                if (std::isfinite(guess) && guess > low && guess <= high) p = guess;
            }
            bool current_is_newton = false;
            bool converged = false;
            for (int iteration = 0; iteration < max_iterations; ++iteration)
            {
                ++record.iterations;
                if (iteration == 0) ++record.initial;
                else if (current_is_newton) ++record.newton;
                else ++record.bisect;
                const auto [f, derivative] = residual(p);
                if (!std::isfinite(f)) throw std::runtime_error("ExactSolver: nonfinite pressure residual");
                if (f > 0) high = p; else low = p;
                if (high - low <= high * tolerance)
                {
                    converged = true;
                    break;
                }
                const double correction = legacy ? 0 : f / derivative;
                if (!legacy && std::isfinite(correction) && std::abs(correction) <= p * tolerance)
                {
                    // Certify a narrow sign-changing bracket: a small Newton step alone
                    // is insufficient near vacuum or when cancellation masks the residual.
                    const double lo = std::max(low, p * (1 - 0.25 * tolerance));
                    const double hi = std::min(high, p * (1 + 0.25 * tolerance));
                    if (lo > 0 && residual(lo).first <= 0 && residual(hi).first >= 0 &&
                        hi - lo <= hi * tolerance)
                    {
                        low = lo; high = hi; converged = true; break;
                    }
                }
                const double candidate = p - correction;
                // Force a bisection periodically to bound stagnation for difficult states.
                if (!legacy && iteration % 8 != 7 && std::isfinite(candidate) &&
                    candidate > low && candidate < high)
                {
                    p = candidate; current_is_newton = true;
                }
                else
                {
                    p = (low + high) / 2;
                    current_is_newton = false;
                }
            }
            if (!converged) throw std::runtime_error("ExactSolver: star pressure did not converge");
            const double pre_star = (high + low) / 2;
            const double f_left = legacy ? wave_condition(pre_star, rho_L, pre_L, gamma_L) : left.evaluate(pre_star).first;
            const double f_right = legacy ? wave_condition(pre_star, rho_R, pre_R, gamma_R) : right.evaluate(pre_star).first;
            const double v_nstar = (v_nL - f_left + v_nR + f_right) / 2;
            record.failed = false;
            return std::pair<double, double>{pre_star, v_nstar};
        }

        /**
         * @brief Sample a vacuum-producing Riemann solution at the stationary interface.
         * @param[in] W_L Original left gas state with positive density and pressure.
         * @param[in] W_R Original right gas state with positive density and pressure.
         * @param[in] gamma_L Left adiabatic index greater than one.
         * @param[in] gamma_R Right adiabatic index greater than one.
         * @param[in] direction Normal axis in [0, hd_core::NDIM).
         * @return Vacuum or a sampled rarefaction state, together with its adiabatic index.
         * @pre The two valid gas states satisfy the vacuum-generation criterion.
         */
        InterfaceState sample_vacuum(const hd_core::Primitives& W_L, const hd_core::Primitives& W_R,double gamma_L, double gamma_R,int direction)
        {
            //deal with vacuum generated in Riemann Problem
            //-----------|~~~~~~~~|   vacuum    |~~~~~~~~|-------------
            const double v_nL = W_L[hd_core::VEL(direction)];
            const double v_nR = W_R[hd_core::VEL(direction)];
            const double c_sL = std::sqrt(
            gamma_L * W_L[hd_core::PRE] / W_L[hd_core::RHO]
            );
            const double c_sR = std::sqrt(
                gamma_R * W_R[hd_core::PRE] / W_R[hd_core::RHO]
            );
            const double S_VL = v_nL + 2.0 * c_sL / (gamma_L - 1.0);//(Left-Vacuum Interace velocity)
            const double S_VR = v_nR - 2.0 * c_sR / (gamma_R - 1.0);
            if(S_VL<=0 && S_VR>=0)
            {
                return {hd_core::Primitives{},gamma_L}; // In the vacuum
            }
            if(S_VL>0.0)
            {
                return {sample_rarefraction(W_L,gamma_L,S_VL,0.0,Side::Left,direction),gamma_L};
            }
            return {sample_rarefraction(W_R,gamma_R,S_VR,0.0,Side::Right,direction),gamma_R};
            //Though there is no star region. In this case, the interface between the tail of rarefraction and vacuum is equivalent to star region, because they are both the boundary of rarefraction(tail).
            //So S_V is v_nstar because it's the speed of boundary, and p_star =0.0, because it can make c_sK_star = 0, the v_rare_tail = S_V....everything fine!
        }
        /**
         * @brief Sample expansion of a gas state into an initially empty neighbor.
         * @param[in] gas Non-vacuum primitive gas state.
         * @param[in] gamma Gas adiabatic index greater than one.
         * @param[in] gas_side Side initially occupied by gas.
         * @param[in] direction Normal axis in [0, hd_core::NDIM).
         * @return Vacuum or a rarefaction state at x/t = 0, with the gas adiabatic index.
         * @pre gas has finite positive density and pressure.
         */
        InterfaceState gas_to_vacuum(const hd_core::Primitives& gas,double gamma, Side gas_side, int direction)
        {
            //To deal with real vacuum cell
            //Basic Idea(take left gas, right vacuum as example)
            //Gas expand to vacuum, producing rarefraction
            //--------------------
            const double v_n = gas[hd_core::VEL(direction)];
            const double c_s = std::sqrt(gamma * gas[hd_core::PRE]/gas[hd_core::RHO]);
            const double S_V = v_n + (gas_side == Side::Left? 1.0 : -1.0) * 2.0 * c_s / (gamma - 1.0);
            if((gas_side == Side::Left && S_V<= 0.0) || (gas_side==Side::Right && S_V>=0.0))
            {
                return {hd_core::Primitives{},gamma};
            }
            return {sample_rarefraction(gas,gamma,S_V,0.0,gas_side,direction),gamma};
            //Though there is no star region. In this case, the interface between the tail of rarefraction and vacuum is equivalent to star region, because they are both the boundary of rarefraction(tail).
            //So S_V is v_nstar because it's the speed of boundary, and p_star =0.0, because it can make c_sK_star = 0, the v_rare_tail = S_V....everything fine!
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
            const hd_core::Primitives& WL,
            const hd_core::Primitives& WR,
            int direction,
            double gamma_L,
            double gamma_R
        ) const
        {
            if (direction < 0 || direction >= hd_core::NDIM)
            {
                throw std::out_of_range(
                    "ExactSolver: invalid direction"
                );
            }

            WL.validate(gamma_L);
            WR.validate(gamma_R);
            const bool vacuum_L = (WL[hd_core::RHO] == 0.0);
            const bool vacuum_R = (WR[hd_core::RHO] == 0.0);
            if (vacuum_L && vacuum_R)
            {
                wave_timestep.face_speed = 0.0;
                return hd_core::Conserved{};
            }
            double min_wave_speed = 0.0;
            double max_wave_speed = 0.0;
            InterfaceState interfaceState;
            if(vacuum_R)
            {
                //left gas to right vacuum
                interfaceState = gas_to_vacuum(WL,gamma_L,Side::Left,direction);
                const double u = WL[hd_core::VEL(direction)];
                const double c = std::sqrt(gamma_L * WL[hd_core::PRE] / WL[hd_core::RHO]);
                min_wave_speed = u - c;
                max_wave_speed = u + 2.0 * c / (gamma_L - 1.0);
            }
            else if(vacuum_L)
            {
                //right gas to left vacuum
                interfaceState = gas_to_vacuum(WR,gamma_R,Side::Right,direction);
                const double u = WR[hd_core::VEL(direction)];
                const double c = std::sqrt(gamma_R * WR[hd_core::PRE] / WR[hd_core::RHO]);
                min_wave_speed = u - 2.0 * c / (gamma_R - 1.0);
                max_wave_speed = u + c;
            }
            else
            {

                auto star = solve_star(WL,WR,gamma_L,gamma_R,tolerance_,max_iterations_,direction,method_);
                if(star.has_value())
                {
                    auto [p_star, vn_star] = *star;

                    // 无论界面采样在哪一侧，都计算左右两侧波速。
                    min_wave_speed = outer_wave_speed(WL, gamma_L, p_star, Side::Left, direction);
                    max_wave_speed = outer_wave_speed(WR, gamma_R, p_star, Side::Right, direction);

                    interfaceState = sample_interface(
                        WL,WR,
                        gamma_L, gamma_R,
                        vn_star, p_star,direction
                    );
                }
                else
                {
                    // 中间生成真空：最外侧边界仍是左右稀疏波头。
                    min_wave_speed = outer_wave_speed(WL, gamma_L, 0.0, Side::Left, direction);
                    max_wave_speed = outer_wave_speed(WR, gamma_R, 0.0, Side::Right, direction);
                    interfaceState = sample_vacuum(
                        WL,WR,
                        gamma_L,gamma_R,
                        direction
                    );
                }
            }
            if (!std::isfinite(min_wave_speed) ||
                !std::isfinite(max_wave_speed) ||
                min_wave_speed > max_wave_speed)
                throw std::runtime_error("ExactSolver: invalid wave-speed bounds");

            // 即使界面通量为零，也把本次最大传播速度交给网格循环。
            wave_timestep.face_speed = std::max(std::abs(min_wave_speed), std::abs(max_wave_speed));
            if (interfaceState.W[hd_core::RHO] == 0.0)
                return hd_core::Conserved{};
            return physical_flux(interfaceState.W, direction, interfaceState.gamma);
        }
}
