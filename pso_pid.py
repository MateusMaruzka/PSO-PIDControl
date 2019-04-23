#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 16:20:09 2019

@author: Mateus Maruzka Roncaglio
"""
import scipy.signal as ss
import pypso_lib as pso
import pypid_control_lib as pid
import pickle


def fObj_pid(x,P,ts,tf,LAMBDA):

    mY, mE, mDU = pid.picontrol(P,ts,tf,x,len(x)) # testar
    mDU = mDU[...,1:-1] - mDU[...,0:-2]
    f = pid.ise(mE) + LAMBDA*pid.ise(mDU) # MULTIOBJETIVO (LQR)
    
    return f


def main():
    
    wmin = 0.1
    wmax = 0.9
    c1 = 2
    c2 = 2
    Ts = 0.1
    Tf = 30
    l = 2
    n_p = 20
    P = ss.TransferFunction([1], [1, 4, 6, 4, 1])
    gb, fitIter = pso.pso(fObj_pid, n_p , 2, _alfa=15, _Wmin=wmin, _Wmax=wmax, _c1 = c1, _c2 = c2, P=P, ts=Ts, tf=Tf, LAMBDA=l)    
    f_name = "dados/pid_pso_n3.pickle"
    
    with open(f_name, "wb") as f:
    
        data = {'Wmin' : wmin,
                'Wmax' : wmax,
                'c1' : c1,
                'c2' : c2,
                'Ts' : Ts,
                'Tf' : Tf,
                'Lambda' : l,
                'Partic' : n_p}
        
        pso_result = {'Gbest': gb,
                      'fitIter' : fitIter } 
        
        pickle.dump(data,f)
        pickle.dump(pso_result, f)

if __name__ == "__main__":

    main()