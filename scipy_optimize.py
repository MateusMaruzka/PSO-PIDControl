# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 21:20:04 2020

@author: maruzka
"""
from scipy import optimize

import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

import pypid_control_lib as pid
import pypso_lib as pso


atraso = 18
P = signal.TransferFunction([1],[2.3335, 1])

    
ts = 0.1
tf = 300

def fObj_pid(x):

    # LAMBDA = 0
    mY, mE, mDU,r,t = pid.picontrol(P,ts,tf,np.array([x]),1, atraso=atraso) # testar
    mDU = mDU[...,1:-1] - mDU[...,0:-2]
    f = np.sum(mE**2) 
    # + 0.2*pid.ise(mDU) # MULTIOBJETIVO (LQR)
    
    return f


def fObj_pid2(x,P,ts,tf,LAMBDA):
    
    #print(P)
    mY, mE, mDU,r,t = pid.picontrol(P[0], ts, tf, x, len(x), atraso = P[1])
    mDU = mDU[...,1:-1] - mDU[...,0:-2]
    f = (1-LAMBDA)*pid.ise(mE) + LAMBDA*pid.ise(mDU) # MULTIOBJETIVO (LQR)
    
    return f

def fObj_pso(x):
    
    wmax, wmin, c1, c2, np = x
    
    atraso = 18;
    P = signal.TransferFunction([1], [2.3335711, 1])
    Ts = 0.1
    Tf = 100
    gb, fitIter = pso.pso(fObj_pid2, int(np) , len(x), var=5, _Wmin=wmin, _Wmax=wmax,
                          _c1 = c1, _c2 = c2, P=[P, atraso], ts=Ts, tf=Tf,
                          LAMBDA=0)  
    
    return fitIter[-1]
    
    
def main():
    
    # sphere = lambda x: np.sum(x**2)
    
    res = optimize.differential_evolution(fObj_pso, ((0.1,0.9), (0.1,0.9), (2.01, 10), (2.01, 10), (10, 100)))
    print(res)
    
    # G = signal.TransferFunction([1],[1,4,6,4,1])

    # mY, mE, mDU,r,t = pid.picontrol(G,ts,tf,np.array([[res[0], res[1]]]),1, atraso = 0) # testar

    # fig, ax = plt.subplots(ncols = 1, nrows = 1)
    
    # ax.plot(mY[0])
    
    

    

if __name__ == "__main__":
    
    main()