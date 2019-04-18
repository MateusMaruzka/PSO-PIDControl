# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 15:59:21 2017

@author: Maruzka
"""

# -*- coding: utf-8 -*-
"""
Exemplo de simulaÃ§Ã£o discreta de uma funÃ§Ã£o de transferÃªncia

Processo de primeira ordem - temperatura ("secador de cabelo")

Y(s)/U(s) = 25/(s + 1)

@author: Prof. Daniel Cavalcanti Jeronymo
"""

import scipy.signal
import numpy as np
import math
import matplotlib.pyplot as plt

# Continuous transfer function of a process
P = scipy.signal.TransferFunction([30303.03], [1, 30303.03])
#P = scipy.signal.TransferFunction([1], [1, 9,23,15])

# Discrete process for simulation
#Ts = 0.01                    # time step


def step(b, a, Ts = 0.1, tf = 1 ):
    
    # Pd = P.to_discrete(Ts)  
    # B = Pd.num                  # zeros
    # A = Pd.den                  # poles
    B = b
    A = a
    nb = len(B) - 1             # number of zeros
    na = len(A) - 1             # number of poles

    # Simulation parameters
    # tf = 0.05
    
    slack = np.amax([na, nb]) + 1 # slack for negative time indexing of arrays
    kend = math.ceil(tf/Ts) + 1   # end of simulation in discrete time
    kmax = kend + slack           # total simulation array size
    
    y = np.zeros(kmax)
   # y.dtype = np.float128
    u = np.ones(kmax)
    #u.dtype = np.float128
    # u[0:slack] = 0
    
    # Simulate
    for k in range(slack, kmax-1):
        #print(len(A[1:]))
        #print(len(y[k-1:k-1-na:-1]))
        y[k] = np.dot(B, np.ones(len(B))) - np.dot(A[1:], y[k-1:k-1-na:-1])

        #y[k] = np.dot(B, u[k-1:k-1-(nb+1):-1]) - np.dot(A[1:], y[k-1:k-1-na:-1])
        
    # ignore empty slack
    y = y[slack:]
    u = u[slack:]
    t = np.arange(0, tf + Ts, Ts)
    
    return y,t

# Plot time response
def main():
    
    Ts = 0.1
    tf = 25
    P = scipy.signal.TransferFunction([1], [1, 4, 6, 4, 1])
    Pd = P.to_discrete(Ts)
    t = np.arange(0, tf + Ts, Ts)
    
    y1,t1 = step(Pd.num, Pd.den, Ts, tf)
    t2,y2 = scipy.signal.dstep(Pd, n=len(t1))
    
    fig, ax = plt.subplots(2, sharex=True)
    
 #   ax[0].plot(t, , 'k--')
    
    ax[0].step(t1, y1)
    ax[0].step(t2, y2[0])
    ax[0].set_ylabel('y(t)')
    
    #ax[1].step(t, u)
    #ax[1].set_ylabel('u(t)')
    
    plt.xlabel('t (s)')
    plt.xlim([0, tf])
    
    plt.show()


if __name__ == "__main__":
    main()



"""

h(z) = y(z)/u(z) = pol(a)/pol(b)

y(z)*pol(b) = pol(a)*u(z)

a(n)*z^(n) + a(n-1)*z^(n-1) +a(n-2)*z^(n-2) ...

b(n)*z^(n) + b(n-1)*z^(n-1) +b(n-2)*z^(n-2) ...

pol(a)*u(z) = [z^(m) + a(1)*z^(m-1) +a(2)*z^(m-2)  + An ...]*[U(z)]

pol(b)*y(z) = [z^(n) + b(1)*z^(n-1) +b(2)*z^(n-2) + ... +Bn ]*[Y(z)]


[z^(n) + b(1)*z^(n-1) +b(2)*z^(n-2) + ... +Bn ]*[Y(z)] = [z^(m) + a(1)*z^(m-1) +a(2)*z^(m-2)  + An ...]*[U(z)]

y[n] + b1*y[n-1] +b2*y[n-2] + ... + bn = u[z] + a1*y[m-1] + a2*z[m-2] + ... + an

y[n] = u[z] + a1*y[m-1] + a2*z[m-2] + ... + an - b1*y[n-1] +b2*y[n-2] + ... + bn

"""

#plt.plot(t, y)
#plt.xlabel('t (s)')
#plt.ylabel('y(t)')
