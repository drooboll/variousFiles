import math

from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt

f0 = 1 * 10**9
# math.sqrt(3) for poly, 20/2 for homework
h = 10
# 20 for poly, 25 for homework
v = 25
# 0 for poly, 1 for homework
d0 = 1
c = 3 * 10**8


def func(tau, t):
    root = math.sqrt((v * t + d0) ** 2 + 4 * h ** 2)

    phi1 = 2 * np.pi * f0 * (v * t + d0) / c
    phi2 = 2 * np.pi * f0 * root / c

    sinc1 = np.sinc(tau - (v * t + d0)/c)
    sinc2 = np.sinc(tau - root / c)
    cos1 = math.cos(phi1)
    sin1 = math.sin(phi1)
    cos2 = math.cos(phi2)
    sin2 = math.sin(phi2)

    # 1 for polycopy - 1 path with reflexion
    # 2 for homework - 2 paths with reflexion
    k = 2
    return (1 / (v * t + d0) * sinc1) ** 2 + k * k * (1 / root * sinc2) ** 2 + 2 * k * ((cos1 * cos2 + sin1 * sin2) / (v * t + d0) / root) * sinc1 * sinc2

listx = []
listy = []

step = math.pow(10 ** 4, 1/1000)
for i in range(0, 1000):
    listx.append(step ** i / 10 ** 3)
    listy.append(integrate.quad(func, -np.inf, np.inf, args=(step ** i / 10 ** 3))[0])
fig, ax = plt.subplots()
ax.plot(listx, listy)

ax.set(xlabel='time (s)', ylabel='Ec', title='2 reflections')
ax.grid()

ax.set_xscale("log", base=10)
ax.set_yscale("log", base=10)
plt.show()

