#pragma once
#include "exact.h"
#include "core.h"
#include <cmath>
#include <optional>
#include <utility>
namespace hydro
{
    namespace
    {
        enum class Side
        {
            Left,
            Right
        };
        double wave_condition(double p_star, double rho_K, double p_K, double gamma_K);
        std::optional<std::pair<double, double>> solve_star(const hd_core::Primitives& W_L, const hd_core::Primitives& W_R,double gamma_L, double gamma_R,double tolerance,int max_iterations,int direction);
        hd_core::Primitives sample_rarefraction(const hd_core::Primitives& K, double gamma_K, double v_nstar, double pre_star,Side side, int direction);
        hd_core::Primitives sample_shock(const hd_core::Primitives& K, double gamma_K, double v_nstar, double pre_star, Side side, int direction);
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
        
        std::optional<std::pair<double, double>> solve_star(const hd_core::Primitives& W_L, const hd_core::Primitives& W_R,double gamma_L, double gamma_R,double tolerance,int max_iterations,int direction)
        {
            double rho_L = W_L[hd_core::RHO];
            double rho_R = W_R[hd_core::RHO];
            double v_nL = W_L[hd_core::VEL(direction)];
            double v_nR = W_R[hd_core::VEL(direction)];
            double pre_L = W_L[hd_core::PRE];
            double pre_R = W_R[hd_core::PRE];
            double cs_L = std::sqrt(gamma_L * pre_L/rho_L);
            double cs_R = std::sqrt(gamma_R * pre_R/rho_R);
            double v_sL = v_nL + 2 * cs_L / (gamma_L - 1);
            double v_sR = v_nR - 2 * cs_R / (gamma_R - 1);
            auto residual = [&](double p) {
                return wave_condition(p, W_L[hd_core::RHO],W_L[hd_core::PRE],gamma_L) + wave_condition(p, W_R[hd_core::RHO],W_R[hd_core::PRE],gamma_R) + v_nR - v_nL;
            };
            if (residual(0.0) >= 0.0)
            {
                return std::nullopt;
            }
            double pre_low = 0;
            double pre_high = std::max(pre_L,pre_R);
            while (residual(pre_high) < 0) pre_high *= 2;
            int count = 0;
            do
            {
                count++;
                double mid = (pre_low + pre_high) /2;
                if(residual(mid)>0.0)
                {
                    pre_high = mid;
                }
                else
                {
                    pre_low = mid;
                }
                /* code */
            } while ((pre_high - pre_low) > pre_high * tolerance && count<max_iterations);
            double pre_star = (pre_high+pre_low)/2;
            double v_nstar = (v_nL - wave_condition(pre_star,rho_L,pre_L, gamma_L) + v_nR + wave_condition(pre_star,rho_R,pre_R,gamma_R))/2;
            return std::pair<double, double>{pre_star, v_nstar};

        }
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
    ExactSolver::ExactSolver(double tolerance, int max_iterations)
    : tolerance_(tolerance),
      max_iterations_(max_iterations)
    {
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
                return hd_core::Conserved{};
            }
            InterfaceState interfaceState;
            if(vacuum_R)
            {
                interfaceState = gas_to_vacuum();
            }
            else if(vacuum_L)
            {
                interfaceState = gas_to_vacuum();
            }
            else
            {

                auto rho_L = WL[hd_core::RHO];
                auto rho_R = WR[hd_core::RHO];

                auto v_nL = WL[hd_core::VEL(direction)];
                auto v_nR = WR[hd_core::VEL(direction)];

                auto p_L = WL[hd_core::PRE];
                auto p_R = WL[hd_core::PRE];

                
                auto star = solve_star(WL,WR,gamma_L,gamma_R,tolerance_,max_iterations_,direction);
                if(star.has_value())
                {
                    auto [p_star, vn_star] = *star;

                    interfaceState = sample_interface(
                        WL,WR,
                        gamma_L, gamma_R,
                        vn_star, p_star,direction
                    );
                }
                else
                {
                    interfaceState = sample_vacuum(
                        WL,WR,
                        gamma_L,gamma_R,
                        direction
                    );
                }
            }
            if (interfaceState.W[hd_core::RHO] == 0.0)
            {
                return hd_core::Conserved{};
            }

            return physical_flux(interfaceState.W, direction, interfaceState.gamma);
        }
}