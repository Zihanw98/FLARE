#include "hd_quantities.h"
#include <cmath>
#include <stdexcept>

Grid::Grid() : Grid(1.0, 0.0, 1.0)
{
    this->gamma = 5.0/3.0;
}

Grid::Grid(double density, double velocity_x, double pressure)
{
    if (!std::isfinite(density) || density <= 0.0)
        throw std::invalid_argument("Density should be finite positive!");
    if (!std::isfinite(velocity_x))
        throw std::invalid_argument("Velocity X should be finite!");
    if (!std::isfinite(pressure) || pressure <= 0.0)
        throw std::invalid_argument("Pressure should be finite positive!");

    this->rho = density;
    this->vx = velocity_x;
    this->p = pressure;
    this->gamma = 5.0/3.0;
}
Grid::Grid(double density, double velocity_x, double pressure,double gamma)
{
    if (!std::isfinite(density) || density <= 0.0)
        throw std::invalid_argument("Density should be finite positive!");
    if (!std::isfinite(velocity_x))
        throw std::invalid_argument("Velocity X should be finite!");
    if (!std::isfinite(pressure) || pressure <= 0.0)
        throw std::invalid_argument("Pressure should be finite positive!");

    this->rho = density;
    this->vx = velocity_x;
    this->p = pressure;
    this->gamma = gamma;
}
Grid::Grid(const State& prim)
{
    this->rho = prim[0];
    this->vx = prim[1];
    this->p = prim[2];
    this->gamma = 5.0/3.0;
}

Grid Grid::from_Conserved(const State& conserved,double gamma)
{
    const double density = conserved[0];
    const double momentum = conserved[1];
    const double energy = conserved[2];

    if (!std::isfinite(density) || density <= 0.0)
        throw std::invalid_argument("Density should be finite positive!");
    if (!std::isfinite(momentum) || !std::isfinite(energy))
        throw std::invalid_argument("Momentum and Energy density should be finite");

    const double velocity_x = momentum / density;
    const double pressure = (gamma- 1.0)
                          * (energy - 0.5 * momentum * velocity_x);
    return Grid(density, velocity_x, pressure,gamma);
}

State Grid::get_Prim() const
{
    return {rho, vx, p};
}

State Grid::get_Conserved() const
{
    const double E = p / (gamma - 1.0) + 0.5 * rho * vx * vx;
    return {rho, rho * vx, E};
}

State Grid::get_Flux() const
{
    const double E = get_Conserved()[2];
    return {rho * vx, rho * vx * vx + p, (E + p) * vx};
    //Grid Flux, not interface value
}
