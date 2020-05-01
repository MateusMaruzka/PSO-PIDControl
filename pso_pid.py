#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 16:20:09 2019

@author: Mateus Maruzka Roncaglio
"""
# import numpy as np
import scipy.signal as ss
import pypso_lib as pso
import pypid_control_lib as pid
import pickle

from functools import wraps

func = pid.itae
@wraps(func)
def fObj_pid(x,P,ts,tf,LAMBDA):
    
    #print(P)
    mY, mE, mDU,r,t = pid.picontrol(P[0], ts, tf, x, len(x), atraso = P[1])
    mDU = mDU[...,1:-1] - mDU[...,0:-2]
    f = (1-LAMBDA)*func(mE) + LAMBDA*pid.ise(mDU) # MULTIOBJETIVO (LQR)
    
    return f

def main():
    
    # np.random.seed(0)
    
    wmin = 0.1
    wmax = 0.9
    c1 = 2.05
    c2 = 2.05
    Ts = 0.1
    Tf = 100
    l = 0.9995
    n_p = 100
    
    
    atraso = 18;
    P = ss.TransferFunction([1], [2.3335711, 1])
    # [18.         2.3335711]
    # P = ss.TransferFunction([1],[1, 4, 6, 4, 1])

    gb, fitIter = pso.pso(fObj_pid, n_p , 2, var=5, _Wmin=wmin, _Wmax=wmax,
                          _c1 = c1, _c2 = c2, P=[P, atraso], ts=Ts, tf=Tf,
                          LAMBDA=l)    


 
    metodo = "{:s}_l={:2.4f}".format(fObj_pid.__name__.upper(), l)

    f_name = "dados/"+metodo+".pickle"
    
    print(gb, fitIter[-1])
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
    
    
"""
 fun: 26.35480948668638
 message: 'Optimization terminated successfully.'
    nfev: 264
     nit: 1
 success: True
       x: array([ 0.11625663,  0.12544062,  2.38722146,  5.90464252, 98.53280363])

"""