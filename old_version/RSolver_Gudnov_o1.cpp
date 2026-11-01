#pragma once
#include "RiemannSolver.h"
#include <iostream>
#include <math.h>

/**
 * @brief Get Flux in the interface
 * 
 * @param grid_L the left grid
 * @param grid_R the right grid
 * @
 */
State Godunovo1_GetInterfaceFlux(const Grid& grid_L, const Grid& grid_R)
{
    auto [pre_star,vx_star] = Godunovo1_SolveStar(grid_L.gamma,grid_R.gamma,grid_L.get_Prim(),grid_R.get_Prim());
    State prim_L,prim_R;
    prim_L = grid_L.get_Prim();
    prim_R = grid_R.get_Prim();
    double rho_L = prim_L[0],rho_R = prim_R[0];
    double vx_L = prim_L[1],vx_R = prim_R[1];
    double pre_L = prim_L[2], pre_R = prim_R[2];
    double rho_f,vx_f,pre_f;
    if(vx_star > 0)     //Grid interface is left of the contact interface--->shock or rarefraction
    {
        const double g = grid_L.gamma;
        double beta = (g - 1) / (g + 1);
        double c_L = std::sqrt(g * pre_L/rho_L); 
        double comp_ratio = pre_star / pre_L;
        double rho_star;
        if(comp_ratio>1)
        {
            //left shock
            double shock_speed = vx_L -c_L*std::sqrt((g+1)/(2*g)*comp_ratio+(g-1)/(2*g));
            if(shock_speed>0)
            {
                //Grid interface is left of the shock front...So interface is at the same status as before
                rho_f = prim_L[0];
                vx_f = prim_L[1];
                pre_f = prim_L[2];
            }
            else
            {
                //star region
                rho_star = rho_L * (comp_ratio + beta) / ( beta * comp_ratio + 1);
                rho_f = rho_star;
                vx_f = vx_star;
                pre_f = pre_star;
            }
        }
        else
        {
            //left rarefraction
            double c_L_star = c_L * std::pow(pre_star/pre_L,(g-1)/(2*g));
            if(vx_L - c_L >= 0)
            {
                // interface at original left region
                rho_f = prim_L[0];
                vx_f = prim_L[1];
                pre_f = prim_L[2];
                
            }
            else
            {
                if(vx_star - c_L_star<=0)
                {
                    rho_star = rho_L * std::pow(pre_star / pre_L, 1.0 / g);
                    rho_f = rho_star;
                    vx_f = vx_star;
                    pre_f = pre_star;
                    //interface at the left star region
                }
                else
                {
                    vx_f = 1/(g+1)*((g-1)*vx_L + 2*c_L);
                    double c_f = vx_f;
                    rho_f = rho_L * std::pow((c_f/c_L),2 / (g - 1));
                    pre_f = pre_L * std::pow((c_f/c_L),2 * g / (g-1));
                    //interface inside the rarefraction wave
                }
            }

        }
    }
    else
    {
        const double g = grid_R.gamma;
        double beta = (g - 1) / (g + 1);
        double c_R = std::sqrt(g * pre_R/rho_R); 
        double comp_ratio = pre_star / pre_R;
        double rho_star;
        if(comp_ratio>1)
        {
            //right shock
            double shock_speed = vx_R + c_R *std::sqrt((g+1)/(2*g)*comp_ratio+(g-1)/(2*g));
            if(shock_speed<0)
            {
                rho_f = prim_R[0];
                vx_f = prim_R[1];
                pre_f = prim_R[2];
            }
            else
            {
                //star region
                rho_star = rho_R * (comp_ratio + beta) / ( beta * comp_ratio + 1);
                rho_f = rho_star;
                vx_f = vx_star;
                pre_f = pre_star;
            }
        }
        else
        {
            double c_R_star = c_R * std::pow(pre_star/pre_R,(g-1)/(2*g));
            if(vx_R + c_R <= 0)
            {
                rho_f = prim_R[0];
                vx_f = prim_R[1];
                pre_f = prim_R[2];
            }
            else
            {
                if(vx_star + c_R_star >=0)
                {
                    rho_star = rho_R * std::pow(pre_star / pre_R, 1.0 / g);
                    rho_f = rho_star;
                    vx_f = vx_star;
                    pre_f = pre_star;
                }
                else
                {
                    vx_f = 1/(g+1)*((g-1)*vx_R - 2*c_R);
                    double c_f = -vx_f;
                    rho_f = rho_R * std::pow((c_f/c_R),2 / (g - 1));
                    pre_f = pre_R * std::pow((c_f/c_R),2 * g / (g-1));
                }
            }
        }
    }
    double gamma_f = (vx_star > 0.0) ? grid_L.gamma : grid_R.gamma;
    double E_f = pre_f / (gamma_f - 1.0)
         + 0.5 * rho_f * vx_f * vx_f;
    return {rho_f * vx_f, rho_f * vx_f * vx_f + pre_f, (E_f + pre_f) * vx_f};
}
/**
 * @brief Bi-section search to solve the velocity and pressure in star region between grids
 *
 * @param gamma_L adiabatic index of left region
 * @param gamma_R adiabitic index of right region
 * @param prim_L Tri-tuple of prim state of left region
 * @param prim_R Tri-tuple of prim state of right region
 *
 * @return a pair of star region's velocity and pressure
 */
std::pair<double, double>  Godunovo1_SolveStar(double gamma_L,double gamma_R,const State& prim_L, const State& prim_R)
{
    double rho_L = prim_L[0];
    double vx_L = prim_L[1];
    double p_L = prim_L[2];

    double rho_R = prim_R[0];
    double vx_R = prim_R[1];
    double p_R = prim_R[2];

    auto residual = [&](double p) {
        return Godunovo1_WaveCondition(p, gamma_L, prim_L) + Godunovo1_WaveCondition(p, gamma_R, prim_R) + vx_R - vx_L;
    };
    if (residual(0.0) >= 0)
        throw std::runtime_error("This code can't deal with vacuum solution...");
    double pre_low = 0;
    double pre_high = std::max(prim_L[2],prim_R[2]);
    while (residual(pre_high) < 0) pre_high *= 2;
    int count = 0, maxiter = 100;
    const double relative_tolerance = 1e-8;
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
    } while ((pre_high - pre_low) > pre_high * relative_tolerance && count<maxiter);
    double pre_star = (pre_high+pre_low)/2;
    double vx_star = (vx_L - Godunovo1_WaveCondition(pre_star,gamma_L,State{rho_L,vx_L,p_L}) + vx_R + Godunovo1_WaveCondition(pre_star,gamma_R,State{rho_R, vx_R, p_R}))/2;
    return {pre_star, vx_star};
    
}

/**
 * @brief Get the connecting (star) region's velocity delta function f_K(p*). 
 * 
 * @param p_star The pressure of the star region. Should be a same value connecting both sides.
 * @param K The tri-tuple of prim state of contact region(either left or right). Sequence: Density, X Velocity, Pressure
 * 
 * 
 */
double Godunovo1_WaveCondition(double p_star, double gamma_K, const State& K)
{
    double ret = 0.0;
    double rho_K = K[0];
    double p_K = K[2];
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