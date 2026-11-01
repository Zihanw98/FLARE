import numpy as np
import matplotlib.pyplot as plt

# 每一列分别为 x、密度、速度、压力
data = np.loadtxt("sod_cpp.csv", delimiter=",", skiprows=1)
x, rho, u, p = data.T

# 将一维密度复制成二维数组：(60, N)
tube = np.tile(rho, (60, 1))

# x 是单元中心，绘图范围需要延伸到单元边界
dx = x[1] - x[0]
xmin = x[0] - dx / 2
xmax = x[-1] + dx / 2
height = 0.12

fig, ax = plt.subplots(figsize=(10, 2.5), layout="constrained")

im = ax.imshow(
    tube,
    origin="lower",
    extent=[xmin, xmax, 0, height],
    aspect="equal",
    cmap="viridis",
    interpolation="nearest",
)

ax.set_xlabel("x")
ax.set_yticks([])
ax.set_title("Sod shock tube: density at t = 0.2")
fig.colorbar(im, ax=ax, label="Density", pad=0.02)

fig.savefig("sod_tube.png", dpi=200)
plt.show()
