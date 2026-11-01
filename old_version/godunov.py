"""一维 Euler 方程：一阶 Godunov + 精确 Riemann 求解器。
只依赖 Python 标准库。运行：python3 godunov.py
输出当前目录下的 sod_python.csv：x, rho, u, p。
教学版本：理想气体，无源项，不处理真空。
"""
import csv
from math import sqrt

GAMMA = 1.4


def conserved(W):
    """W=(rho,u,p) -> U=(rho,rho*u,E)。"""
    rho, u, p = W
    return [rho, rho*u, p/(GAMMA-1) + 0.5*rho*u*u]


def primitive(U):
    rho, m, E = U
    if rho <= 0:
        raise ValueError("密度必须为正")
    u = m/rho
    p = (GAMMA-1)*(E - 0.5*m*u)
    if p <= 0:
        raise ValueError("压力必须为正")
    return rho, u, p


def physical_flux(W):
    rho, u, p = W
    E = p/(GAMMA-1) + 0.5*rho*u*u
    return [rho*u, rho*u*u+p, (E+p)*u]


def wave_function(p, W):
    """压力方程中的 f_K(p)，由激波/稀疏波关系得到。"""
    rho, u, p0 = W
    c = sqrt(GAMMA*p0/rho)
    if p > p0:  # 激波
        A = 2/((GAMMA+1)*rho)
        B = (GAMMA-1)/(GAMMA+1)*p0
        return (p-p0)*sqrt(A/(p+B))
    # 稀疏波
    return 2*c/(GAMMA-1)*((p/p0)**((GAMMA-1)/(2*GAMMA))-1)


def star_pressure_velocity(L, R):
    """用二分法求 f_L(p*)+f_R(p*)+u_R-u_L=0。"""
    def residual(p):
        return wave_function(p, L) + wave_function(p, R) + R[1]-L[1]

    if residual(0.0) >= 0:
        raise ValueError("这个教学版本不处理真空 Riemann 解")
    lo, hi = 0.0, max(L[2], R[2])
    while residual(hi) < 0:
        hi *= 2
    for _ in range(60):
        mid = 0.5*(lo+hi)
        if residual(mid) > 0:
            hi = mid
        else:
            lo = mid
    ps = 0.5*(lo+hi)
    us = 0.5*(L[1]+R[1] + wave_function(ps, R)-wave_function(ps, L))
    return ps, us


def interface_state(L, R, ps, us):
    """精确 Riemann 解在 xi=(x-x_face)/(t-t_n)=0 处的状态。"""
    g = GAMMA
    beta = (g-1)/(g+1)
    if us >= 0:  # 界面位于接触间断左侧
        rho, u, p = L
        c = sqrt(g*p/rho)
        ratio = ps/p
        if ps > p:  # 左激波
            speed = u-c*sqrt((g+1)/(2*g)*ratio+(g-1)/(2*g))
            if speed >= 0:
                return L
            return rho*(ratio+beta)/(beta*ratio+1), us, ps
        # 左稀疏波：依次判断原始区、星区、扇形内部
        cs = c*ratio**((g-1)/(2*g))
        if u-c >= 0:
            return L
        if us-cs <= 0:
            return rho*ratio**(1/g), us, ps
        uf = 2/(g+1)*(c+0.5*(g-1)*u)
        cf = uf  # xi=0 处 u-c=0
    else:  # 界面位于接触间断右侧
        rho, u, p = R
        c = sqrt(g*p/rho)
        ratio = ps/p
        if ps > p:  # 右激波
            speed = u+c*sqrt((g+1)/(2*g)*ratio+(g-1)/(2*g))
            if speed <= 0:
                return R
            return rho*(ratio+beta)/(beta*ratio+1), us, ps
        cs = c*ratio**((g-1)/(2*g))
        if u+c <= 0:
            return R
        if us+cs >= 0:
            return rho*ratio**(1/g), us, ps
        uf = 2/(g+1)*(-c+0.5*(g-1)*u)
        cf = -uf  # xi=0 处 u+c=0
    return rho*(cf/c)**(2/(g-1)), uf, p*(cf/c)**(2*g/(g-1))


def godunov_flux(UL, UR):
    """UL, UR -> 精确界面状态 -> F。"""
    L, R = primitive(UL), primitive(UR)
    ps, us = star_pressure_velocity(L, R)
    Wf = interface_state(L, R, ps, us)
    return physical_flux(Wf)


def solve(n=200, t_end=0.2, cfl=0.4):
    dx = 1.0/n
    U = [conserved((1.0, 0.0, 1.0) if (i+0.5)*dx < 0.5
                   else (0.125, 0.0, 0.1)) for i in range(n)]
    t = 0.0
    while t < t_end:
        W = [primitive(q) for q in U]
        amax = max(abs(u)+sqrt(GAMMA*p/rho) for rho, u, p in W)
        dt = min(cfl*dx/amax, t_end-t)

        # 零梯度边界：两端各复制一个单元作为 ghost cell。
        ext = [U[0]] + U + [U[-1]]
        # 一阶重构：界面左右状态就是相邻单元的守恒量。
        F = [godunov_flux(ext[j], ext[j+1]) for j in range(n+1)]
        # F[i] 是单元 i 左界面，F[i+1] 是其右界面。
        # 必须同时更新：所有通量都来自旧时刻。
        U = [[U[i][k]-dt/dx*(F[i+1][k]-F[i][k])
              for k in range(3)] for i in range(n)]
        t += dt
    return U


if __name__ == "__main__":
    L = conserved((1.0, 0.0, 1.0))
    R = conserved((0.125, 0.0, 0.1))
    print("初始间断处的 Godunov 通量:", godunov_flux(L, R))
    U = solve()
    with open("sod_python.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["x", "rho", "u", "p"])
        for i, q in enumerate(U):
            writer.writerow([(i+0.5)/len(U), *primitive(q)])
    print("已输出 sod_python.csv，t=0.2")
