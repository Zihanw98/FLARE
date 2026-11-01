// 一维 Euler：一阶 Godunov + 精确 Riemann 求解器。
// 编译：g++ -std=c++17 -O2 godunov.cpp -o godunov
// 运行：./godunov，输出当前目录下的 sod_cpp.csv。
// 教学版本：理想气体，无源项，不处理真空。
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

constexpr double GAMMA = 1.4;
using State = std::array<double, 3>; // U=(rho, rho*u, E)
struct Primitive { double rho, u, p; };

State conserved(Primitive w) {
    return {w.rho, w.rho*w.u, w.p/(GAMMA-1)+0.5*w.rho*w.u*w.u};
}

Primitive primitive(State U) {
    if (U[0] <= 0) throw std::runtime_error("密度必须为正");
    double u = U[1]/U[0];
    double p = (GAMMA-1)*(U[2]-0.5*U[1]*u);
    if (p <= 0) throw std::runtime_error("压力必须为正");
    return {U[0], u, p};
}

State physical_flux(Primitive w) {
    double E = w.p/(GAMMA-1)+0.5*w.rho*w.u*w.u;
    return {w.rho*w.u, w.rho*w.u*w.u+w.p, (E+w.p)*w.u};
}

// 压力方程中的 f_K(p)。
double wave_function(double p, Primitive w) {
    double c = std::sqrt(GAMMA*w.p/w.rho);
    if (p > w.p) { // 激波
        double A = 2/((GAMMA+1)*w.rho);
        double B = (GAMMA-1)/(GAMMA+1)*w.p;
        return (p-w.p)*std::sqrt(A/(p+B));
    }
    // 稀疏波
    return 2*c/(GAMMA-1)*(std::pow(p/w.p, (GAMMA-1)/(2*GAMMA))-1);
}

std::pair<double, double> star_pressure_velocity(Primitive L, Primitive R) {
    auto residual = [&](double p) {
        return wave_function(p, L)+wave_function(p, R)+R.u-L.u;
    };
    if (residual(0.0) >= 0)
        throw std::runtime_error("这个教学版本不处理真空 Riemann 解");
    double lo = 0, hi = std::max(L.p, R.p);
    while (residual(hi) < 0) hi *= 2;
    // 二分法求星区压力。
    for (int j = 0; j < 60; ++j) {
        double mid = 0.5*(lo+hi);
        if (residual(mid) > 0) hi = mid;
        else lo = mid;
    }
    double ps = 0.5*(lo+hi);
    double us = 0.5*(L.u+R.u+wave_function(ps, R)-wave_function(ps, L));
    return {ps, us};
}

// 采样精确 Riemann 解在 xi=(x-x_face)/(t-t_n)=0 处的状态。
Primitive interface_state(Primitive L, Primitive R, double ps, double us) {
    const double g = GAMMA, beta = (g-1)/(g+1);
    Primitive w;
    double c, uf, cf;
    if (us >= 0) { // 接触间断左侧
        w = L;
        c = std::sqrt(g*w.p/w.rho);
        double ratio = ps/w.p;
        if (ps > w.p) { // 左激波
            double speed = w.u-c*std::sqrt((g+1)/(2*g)*ratio+(g-1)/(2*g));
            if (speed >= 0) return L;
            return {w.rho*(ratio+beta)/(beta*ratio+1), us, ps};
        }
        // 左稀疏波：原始区、星区、扇形内部。
        double cs = c*std::pow(ratio, (g-1)/(2*g));
        if (w.u-c >= 0) return L;
        if (us-cs <= 0) return {w.rho*std::pow(ratio, 1/g), us, ps};
        uf = 2/(g+1)*(c+0.5*(g-1)*w.u);
        cf = uf;
    } else { // 接触间断右侧
        w = R;
        c = std::sqrt(g*w.p/w.rho);
        double ratio = ps/w.p;
        if (ps > w.p) { // 右激波
            double speed = w.u+c*std::sqrt((g+1)/(2*g)*ratio+(g-1)/(2*g));
            if (speed <= 0) return R;
            return {w.rho*(ratio+beta)/(beta*ratio+1), us, ps};
        }
        double cs = c*std::pow(ratio, (g-1)/(2*g));
        if (w.u+c <= 0) return R;
        if (us+cs >= 0) return {w.rho*std::pow(ratio, 1/g), us, ps};
        uf = 2/(g+1)*(-c+0.5*(g-1)*w.u);
        cf = -uf;
    }
    return {w.rho*std::pow(cf/c, 2/(g-1)), uf,
            w.p*std::pow(cf/c, 2*g/(g-1))};
}

State godunov_flux(State UL, State UR) {
    Primitive L = primitive(UL), R = primitive(UR);
    auto [ps, us] = star_pressure_velocity(L, R);
    return physical_flux(interface_state(L, R, ps, us));
}

int main() {
    const int N = 200;
    const double dx = 1.0/N, cfl = 0.4, t_end = 0.2;
    std::vector<State> U(N), next(N), F(N+1);
    for (int i = 0; i < N; ++i)
        U[i] = conserved((i+0.5)*dx < 0.5 ? Primitive{1,0,1}
                                                 : Primitive{0.125,0,0.1});
    State example = godunov_flux(conserved({1,0,1}), conserved({0.125,0,0.1}));
    std::cout << std::setprecision(10) << "Initial Godunov Flux: ";
    for (double f : example) std::cout << f << ' ';
    std::cout << '\n';

    double t = 0;
    while (t < t_end) {
        double amax = 0;
        for (State q : U) {
            Primitive w = primitive(q);
            amax = std::max(amax, std::abs(w.u)+std::sqrt(GAMMA*w.p/w.rho));
        }
        double dt = std::min(cfl*dx/amax, t_end-t);
        // Zero Gradient Boundary(copy to ghost cell)
        F[0] = godunov_flux(U[0], U[0]);
        for (int j = 1; j < N; ++j)
            F[j] = godunov_flux(U[j-1], U[j]); 
        F[N] = godunov_flux(U[N-1], U[N-1]);
        // Update Flux
        for (int i = 0; i < N; ++i)
            for (int k = 0; k < 3; ++k)
                next[i][k] = U[i][k]-dt/dx*(F[i+1][k]-F[i][k]);
        U.swap(next);
        t += dt;
        
    }

    std::ofstream out("sod_cpp.csv");
    if (!out) throw std::runtime_error("无法创建 sod_cpp.csv");
    out << std::setprecision(17) << "x,rho,u,p\n";
    for (int i = 0; i < N; ++i) {
        Primitive w = primitive(U[i]);
        out << (i+0.5)*dx << ',' << w.rho << ',' << w.u << ',' << w.p << '\n';
    }
    std::cout << "Output to sod_cpp.csv，t=0.2\n";
}
