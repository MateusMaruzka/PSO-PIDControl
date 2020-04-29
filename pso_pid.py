#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 16:20:09 2019

@author: Mateus Maruzka Roncaglio
"""
import numpy as np
import scipy.signal as ss
import pypso_lib as pso
import pypid_control_lib as pid
import pickle


def fObj_pid(x,P,ts,tf,LAMBDA):

    #print(P)
    mY, mE, mDU,r,t = pid.picontrol(P[0],ts,tf,x,len(x), atraso = P[1]) # testar
    # mDU = mDU[...,1:-1] - mDU[...,0:-2]
    f = (1-LAMBDA)*pid.ise(mE) + LAMBDA*pid.ise(mDU) # MULTIOBJETIVO (LQR)
    
    return f




def main():
    
    np.random.seed(0)
    
    wmin = 0.1
    wmax = 0.2
    c1 = 2.05
    c2 = 2.05
    Ts = 0.1
    Tf = 200
    l = 0
    n_p = 10
    
    
    atraso = 18;
    P = ss.TransferFunction([1], [2.3335711, 1])
    # [18.         2.3335711]
    # P = ss.TransferFunction([1],[1, 4, 6, 4, 1])

    gb, fitIter = pso.pso(fObj_pid, n_p , 2, var=2, _Wmin=wmin, _Wmax=wmax,
                          _c1 = c1, _c2 = c2, P=[P, atraso], ts=Ts, tf=Tf,
                          LAMBDA=l)    

    metodo = "ISE_l={:2.4f}".format(l)
    f_name = "dados/"+metodo+".pickle"
    
    with open(f_name, "wb") as f:
    
        data = {'Wmin' : wmin,
                'Wmax' : wmax,
                'c1' : c1,
                'c2' : c2,
                'Ts' : Ts,
                'Tf' : Tf,
                'Lambda' : l,
                'Partic' : n_p}
        
        pso_result = {'Metodo' : metodo,
                      'Process': [P, atraso], 
                      'Gbest': gb,
                      'Params' : data,
                      'fitIter' : fitIter } 
        

        pickle.dump(pso_result, f)

if __name__ == "__main__":

    main()
    print("done")