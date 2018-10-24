import numpy as np
from matplotlib import mlab
import mod1 as mod
import matplotlib.pyplot as plt

Rs = 100
beta = 7.6349*10**(-7)
Vs = [-1, 0, 0]
interv = [-2.55, 2.57]
mu = 1.254693
n = int(1e5)

    arr = np.zeros_like(rx)
    N = len(rx)
    for i in range(N):
        print(rx[i])
        r = [rx[i], 0, 0]
        arr[i] = mod.conc(r, interv, Vs, mu, beta, Rs, 1e-17)
    return arr

x = np.linspace(0.5, 100, 50)
y1 = con(x)
y2 = con(-x)
fig, ax = plt.subplots()                        
ax.plot(x, y1, color="lime", label="upwind")
ax.plot(x, y2, color="blue", label="downwind")
ax.set_xlabel("r")                              
ax.set_ylabel("n(r)") 
ax.legend() 
plt.show() 


