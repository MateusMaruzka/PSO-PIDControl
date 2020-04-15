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
    mY, mE, mDU,r,t = pid.picontrol(P[0],ts,tf,x,len(x)) # testar
    mDU = mDU[...,1:-1] - mDU[...,0:-2]
    f = pid.ise(mE) + LAMBDA*pid.ise(mDU) # MULTIOBJETIVO (LQR)
    
    return f




def main():
    
    wmin = 0.1
    wmax = 0.9
    c1 = 2.05
    c2 = 2.05
    Ts = 0.1
    Tf = 90
    l = 0
    n_p = 100
    
    
    atraso = 17;
    P = ss.TransferFunction([1], [2.335, 1])
    
    gb, fitIter = pso.pso(fObj_pid, n_p , 2, var=30, _Wmin=wmin, _Wmax=wmax,
                          _c1 = c1, _c2 = c2, P=[P, atraso], ts=Ts, tf=Tf,
                          LAMBDA=l)    

    metodo = "ISE"
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