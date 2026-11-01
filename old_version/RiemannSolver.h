#pragma once
#include <utility>
#include "hd_quantities.h"
std::pair<double, double>  Godunovo1_SolveStar(double gamma_L,double gamma_R,const State& prim_L, const State& prim_R);
double Godunovo1_WaveCondition(double p_star, double gamma_K, const State& K);
State Godunovo1_GetInterfaceFlux(const Grid&, const Grid&);
