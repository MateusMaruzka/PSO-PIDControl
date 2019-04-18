# -*- coding: utf-8 -*-
"""
Exemplo de controle PI (proporcional integral)

Processo de primeira ordem - temperatura ("secador de cabelo")
Y(s)/U(s) = 25/(s + 1)

@author: Prof. Daniel Cavalcanti Jeronymo
"""

import scipy.signal
import numpy as np
import math
import matplotlib.pyplot as plt

# Continuous transfer function of a process
P = scipy.signal.TransferFunction([1], [1, 4, 6, 4, 1])

# Discrete process for simulation
Ts = 0.5                    # time step
Pd = P.to_discrete(Ts)

B = Pd.num                  # zeros
A = Pd.den                  # poles
nb = len(B) - 1             # number of zeros
na = len(A) - 1             # number of poles


# Simulation parameters
tf = 30

slack = np.amax([na, nb]) + 1 # slack for negative time indexing of arrays
kend = math.ceil(tf/Ts) + 1   # end of simulation in discrete time
kmax = kend + slack           # total simulation array size

y = np.zeros(kmax)
u = np.zeros(kmax)
e = np.zeros(kmax)
r = 1*np.ones(kmax)


k = 1
tau = 2.5128
ts = 1.8463
theta = ts
t = tau

# ZIEGLER NICHOLS
kp = 0.9*tau/ts
ti = ts/0.3

# SIMC
kp = (2*t + theta)/(3*theta*k)
ti = min(t+theta/2, 8*theta)

# COHEN COON
kp = (t/theta)*(0.9 + theta/(12*t))/k
ti = theta*(30+3/t)/(13+8/t)

# Ver outros como:
# Chien Hrones Reswick
# Astrom Hagglund
# AMIGO

ki = kp/ti

# Simulate
for k in range(slack, kmax):
    y[k] = np.dot(B, u[k-1:k-1-(nb+1):-1]) - np.dot(A[1:], y[k-1:k-1-na:-1])
    
    # error
    e[k] = r[k] - y[k]
    
    # PI control discretized by backwards differences, Kp=0.2, Ki = 0.5
    du = kp*(e[k] - e[k-1])  + ki*e[k]*Ts
    u[k] = u[k-1] + du
    
# ignore empty slack
y = y[slack:]
u = u[slack:]
e = e[slack:]
r = r[slack:]

# Plot time response


t = np.arange(0, tf + Ts, Ts)

t2,y2 = scipy.signal.dstep((Pd.num, Pd.den, Ts), t=t)
fig, ax = plt.subplots(2, sharex=True)
ax[0].plot(t, r, 'k--')
ax[0].step(t, y)
ax[0].step(t2, y2[0])

ax[0].set_ylabel('y(t)')
ax[1].step(t, u)
ax[1].set_ylabel('u(t)')
plt.xlabel('t (s)')
plt.xlim([0, tf])
plt.show()


#f __main__ == "__main__":
    
    
    