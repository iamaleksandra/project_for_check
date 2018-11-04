# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 14:03:55 2018

@author: aleksandra
"""
import mod1 as mod
import numpy as np
import matplotlib.pyplot as plt


def makeData(init, final, step):
    x = np.arange(init, final, step)
    y = np.arange(init, final, step)
    xgrid, ygrid = np.meshgrid(x, y)
    return xgrid, ygrid


Rs = 100
beta = 1e-7
Vs = np.array([0, 0, -1])
Vin, dV, Vfin = -2.45, 0.06, 2.47
N = int((Vfin-Vin)/dV)
Vfin = Vin + N*dV
mu = 1.2

r = np.array([0, 0, 1])

arr = np.zeros((N, N))
for i in range(N):
    u = Vin + i*dV
    for j in range(N):
        w = Vin + j*dV
        v = np.array([u, 0, w])
        arr[i][j] = mod.distrfunc2(r, v, Vs, mu, beta, Rs)

xx, yy = makeData(Vin, Vfin, dV)
zz = arr
plt.figure()
cp = plt.contourf(xx, yy, zz, cmap='jet')
plt.colorbar(cp)
plt.title(r'Distribution, 1AU upwind, $\mu$ = {mu}'.format(mu=mu))
plt.xlabel(r'$V_z/V_s$')
plt.ylabel(r'$V_x/V_s$')
plt.show()
