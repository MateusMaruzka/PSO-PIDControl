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


def esfera(x):
    return np.sum((x)**2, axis=1)

def main():
    
    wmin = 0.1
    wmax = 0.9
    c1 = 2
    c2 = 4
    Ts = 0.1
    Tf = 30
    l = 0
    n_p = 50
    P = ss.TransferFunction([1], [1, 4, 6, 4, 1])
    #P = ss.TransferFunction([1], [12, 4])

    gb, fitIter = pso.pso(fObj_pid, n_p , 2, _alfa=10, _Wmin=wmin, _Wmax=wmax, _c1 = c1, _c2 = c2, P=P, ts=Ts, tf=Tf, LAMBDA=l)    
    #gb, fitIter = pso.pso(esfera, n_p, 3)
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
                      'Process': P, 
                      'Gbest': gb,
                      'Params' : data,
                      'fitIter' : fitIter } 
        
        #pickle.dump(data,f)
        pickle.dump(pso_result, f)

if __name__ == "__main__":

    main()