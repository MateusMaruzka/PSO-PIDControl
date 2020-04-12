#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 22:37:02 2019

@author: maruzka
"""

import scipy.signal as ss
import pypso_lib as pso
import pypid_control_lib as pid


import matplotlib.pyplot as plt
import numpy as np





def fObj_pid(x,P,ts,tf,LAMBDA):

    #print(P)
    mY, mE, mDU,r,t = pid.picontrol(P[0],ts,tf,x,len(x)) # testar
    mDU = mDU[...,1:-1] - mDU[...,0:-2]
    f = pid.ise(mE) + LAMBDA*pid.ise(mDU) # MULTIOBJETIVO (LQR)
    
    return f


def esfera(x):
    return np.sum((x)**2, axis=1)

def main():
    
    wmin = 0.5
    wmax = 0.9
    c1 = 2
    c2 = 2
    Ts = 0.1
    Tf = 10
    l = 50
    n_p = 20
    P = ss.TransferFunction([0.5], [1, 1.7])
    atraso = 15;
    # P = ss.TransferFunction([1], [2.19216486, 1])
    
    gb, fitIter = pso.pso(fObj_pid, n_p , 2, _alfa=10, _Wmin=wmin, _Wmax=wmax, _c1 = c1, _c2 = c2, P=[P, atraso], ts=Ts, tf=Tf, LAMBDA=l)    
    
    a = [pso.pso(fObj_pid, i*10 + 10, 2, _alfa=10, _Wmin=wmin, _Wmax=wmax, _c1 = c1, _c2 = c2, P=[P, atraso], ts=Ts, tf=Tf, LAMBDA=l) for i in range(10)] 
    
    
    X = [a[i][0][0] for i in range(len(a))]
    Y = [a[i][0][1] for i in range(len(a))]
    Z = [a[i][1][-1] for i in range(len(a))]

    print(X)
    
    print(Y)
    
    print(Z)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
           
    ax.scatter(X, Y, Z, linewidth=0.6, antialiased=True)
    
    # Customize the z axis.
    
    # Add labels to axes
    ax.set_ylabel("Y")
    ax.set_ylim([min(Y),max(Y)])
    
    ax.set_xlabel("X")
    ax.set_xlim([min(X),max(X)])
    
    ax.set_zlabel("Z")
    ax.set_zlim([min(Z),max(Z)])

    
    
    plt.show()
    
    
    #gb, fitIter = pso.pso(esfera, n_p, 3)
#    metodo = "ISE"
#    f_name = "dados/"+metodo+".pickle"
#    
#    with open(f_name, "wb") as f:
#    
#        data = {'Wmin' : wmin,
#                'Wmax' : wmax,
#                'c1' : c1,
#                'c2' : c2,
#                'Ts' : Ts,
#                'Tf' : Tf,
#                'Lambda' : l,
#                'Partic' : n_p}
#        
#        pso_result = {'Metodo' : metodo,
#                      'Process': [P, atraso], 
#                      'Gbest': gb,
#                      'Params' : data,
#                      'fitIter' : fitIter } 
#        
#        #pickle.dump(data,f)
#        pickle.dump(pso_result, f)

if __name__ == "__main__":

    main()