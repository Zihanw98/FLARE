#pragma once
#include <array>

using State = std::array<double, 3>;

class Grid
{
public:
    State get_Prim() const;       // (rho, vx, p)
    State get_Flux() const;       // (rho*vx, rho*vx^2+p, (E+p)*vx)
    State get_Conserved() const;  // (rho, rho*vx, E)

    Grid();
    Grid(double rho, double v_x, double p);
    Grid(double rho, double v_x, double p,double gamma);
    explicit Grid(const State& prim);
    static Grid from_Conserved(const State& conserved, double gamma);
    double gamma = 5.0/3.0;
    double T = 0.0;

private:
    double rho;
    double vx;
    double p;
};
