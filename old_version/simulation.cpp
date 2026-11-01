#include <iostream>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include "hd_quantities.h"
#include "globalvar.h"
#include "initialize.h"
#include "RiemannSolver.h"
#include "simulation.h"
#include "output.h"
void simulation_main()
{
    int step = 0;
    double tcurrent = 0;
    const double cfl = 0.4;
    double amax = 0;
    ApplyBoundaryConditions(); // Check layout before allocating the flux array.
    if (!std::isfinite(dxa) || dxa <= 0.0 ||
        !std::isfinite(tmax) || tmax < 0.0)
        throw std::invalid_argument("Invalid cell width or end time");

    // nxa active cells have nxa+1 faces. Reuse storage each time step.
    std::vector<State> flux(nxa + 1);
    Output(tcurrent, step, tcurrent >= tmax);
    while(tcurrent<tmax)
    {
        ApplyBoundaryConditions();
        if(step==0)
        {
            //std::cout<<"Initial Status:"<<std::endl;
            for(int x=0;x<nx;x++)
            {
                auto Prim = Domain[x].get_Prim();
                //std::cout<<"Grid "<<x<<" Density:"<<Prim[0]<<" Velocity:"<<Prim[1]<<" Pressure:"<<Prim[2]<<std::endl;
            }
            //std::cout<<std::endl;
        }
        amax = 0.0;
        for (int i = nxg; i < nxg + nxa; ++i) {
            const Grid& q = Domain[i];
            State w = q.get_Prim();
            amax = std::max(amax, std::abs(w[1])+std::sqrt(q.gamma*w[2]/w[0]));
        }
        if (!std::isfinite(amax) || amax <= 0.0)
            throw std::runtime_error("Invalid maximum signal speed");
        dt = std::min(cfl*dxa/amax,tmax - tcurrent);
        if (!std::isfinite(dt) || dt <= 0.0 || tcurrent + dt <= tcurrent)
            throw std::runtime_error("Time step cannot advance the simulation");

        // 1. Compute ALL face fluxes using the old-time state.
        // flux[j] lies between Domain[nxg+j-1] and Domain[nxg+j].
        // j=0 and j=nxa use the left and right ghost cells, respectively.
        for (int j = 0; j <= nxa; ++j)
        {
            const int right = nxg + j;
            flux[j] = Godunovo1_GetInterfaceFlux(Domain[right - 1], Domain[right]);
        }

        // 2. Update active cells only: U_new = U_old - dt/dx*(F_right-F_left).
        // All fluxes are already stored, so in-place cell updates are safe.
        for (int j = 0; j < nxa; ++j)
        {
            const int i = nxg + j;
            State U = Domain[i].get_Conserved();
            for (int k = 0; k < 3; ++k)
                U[k] -= (dt / dxa) * (flux[j + 1][k] - flux[j][k]);

            // Keep the cell's gamma (the current Sod test has uniform gamma).
            // A multi-material model will also need composition/gamma transport.
            Domain[i] = Grid::from_Conserved(U, Domain[i].gamma);
        }
        ApplyBoundaryConditions();
        tcurrent += dt;
        step += 1;
        Output(tcurrent, step, tcurrent >= tmax);
        //std::cout<<"Step "<<step<<" t="<<tcurrent<<" dt="<<dt<<std::endl;
    }
    std::cout<<"Finished!"<<std::endl;
}
