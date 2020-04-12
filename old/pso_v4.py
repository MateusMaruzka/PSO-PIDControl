#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 14:17:19 2019

@author: Mateus Maruzka
"""
import math
import scipy.signal
import numpy as np
import pypso_lib as ppl
import matplotlib.pyplot as plt

from pypid_control_lib import ise


def esfera(x):
    return np.sum((x)**2, axis=1)

def picontrol(P, ts, tf, vetor_ganhos, num_controladores):
    
    Pd = P.to_discrete(ts)
    
    B = Pd.num                  # zeros
    A = Pd.den                  # poles
    nb = len(B) - 1             # number of zeros
    na = len(A) - 1             # number of poles
        
    slack = np.amax([na, nb]) + 1 # slack for negative time indexing of arrays
    kend = math.ceil(tf/ts) + 1   # end of simulation in discrete time
    kmax = kend + slack           # total simulation array size
    
    y = np.zeros([kmax, num_controladores])
    u = np.zeros([kmax, num_controladores])
    e = np.zeros([kmax, num_controladores])
    r = 1*np.ones([kmax, num_controladores])
    
    kp = vetor_ganhos[:,0]
    ki = vetor_ganhos[:,1]

    # Simulate
    for k in range(slack, kmax):
        y[k] = np.dot(B, u[k-1:k-1-(nb+1):-1]) - np.dot(A[1:], y[k-1:k-1-na:-1])
        
        # error
        e[k] = r[k] - y[k]
        
        # PI control discretized by backwards differences
        du = kp*(e[k] - e[k-1])  + ki*e[k]*ts
        u[k] = u[k-1] + du 
         
        # SATURACAO
        # u[k] = min(max(0, u[k]), 2)
        # gg = u[k] > 10
        # u[k,gg] = 10
        
    return y.T, e.T, u.T


def func_pid(x,P,ts,tf,LAMBDA):

    mY, mE, mDU = picontrol(P,ts,tf,x,len(x)) # testar
    _mDU = mDU[...,1:-1] - mDU[...,0:-2]
    f = ise(mE) + LAMBDA*ise(_mDU) # MULTIOBJETIVO (LQR)
    
    return f#, mE, mDU, mY

def pso_of_psos(x):
    
    fit = []
    for i in range(len(x)):
        a = x[i]
        # alfa = a[0]
        wmin = a[0] #if a[0] < 1 else 0.1
        #x[i][0] = wmin
        wmax = a[1] #if a[1] < 1 else 0.9
       # x[i][1] = wmax
        c1   = a[2] #if a[2] > 0 else 0.5
       # x[i][2] = c1
        c2   = a[3] #if a[3] > 0 else 0.5
       # x[i][3] = c2
        # n_partic = int(a[5]) + 1
        
        P = scipy.signal.TransferFunction([1], [1, 4, 6, 4, 1])
        Ts = 0.1
        Tf = 20
        fitIter = []
        gb, fitIter = ppl.pso(func_pid, 20 , 2, _alfa=15, _Wmin=wmin, _Wmax=wmax, _c1 = c1, _c2 = c2, P=P, ts=Ts, tf=Tf, LAMBDA=2)    
        fit.append(fitIter.pop())


    return np.array(fit)
    


def main():


    gb, fit = ppl.pso(esfera, 30,3, _Wmin = 0.1)
    #f,e,u,y = func_fitness(np.array([[2.4,2.4],[5.8,5.8]]), P, Ts,Tf, 0)

    print(gb)
    print(fit.pop())
    plt.semilogx(fit)
    
    
if __name__ == "__main__":
    main()