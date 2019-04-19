#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 20:56:59 2019

@author: Mateus Maruzka Roncaglio
"""
import numpy as np
import scipy.signal 
import control.matlab
import control
import matplotlib.pyplot as plt

def main():
    
    Ts = 0.5
    P = scipy.signal.TransferFunction([1], [1, 4, 6, 4, 1])
    """
    Segurador de ordem 0? 
    
    """
    # P = scipy.signal.TransferFunction([1], [1, 2, 1, 0])
    Pd = P.to_discrete(Ts);
    w, H = scipy.signal.dfreqresp(Pd)

    plt.figure()
    plt.plot(H.real,(H.imag), "b")
    plt.plot(H.real, -H.imag, "r")
    plt.scatter(-1, 0, marker = '|')
    
    
    for i in range(len(-H.imag) - 1):
        if (H.imag[i] > 0 and H.imag[i+1] < 0) or (H.imag[i] < 0 and H.imag[i+1] > 0):
            print(i)
    #print(H.imag[2841])
   # print(H.real[2841])
    j = next(i for i in range(0,len(H.imag)-1) if (H.imag[i] > 0 and H.imag[i+1] < 0) or (H.imag[i] < 0 and H.imag[i+1] > 0))
    a = H.real[j]
    b = H.imag[j]
    # print (1/H.real[next(i for i in range(0,len(H.imag)-1) if (H.imag[i] > 0 and H.imag[i+1] < 0) or (H.imag[i] < 0 and H.imag[i+1] > 0))])
    
    plt.scatter(a, b,c='k',marker = 'o')

   
    # sys = control.matlab.tf(1, [1, 4, 6, 4, 1])
    # sys = control.matlab.tf([1], [1, 2, 1, 0])
    # gm, pm, wg, wp = control.matlab.margin(sys)
    # print(control.matlab.margin(sys))

    plt.show()

      



if __name__ == "__main__":
    main()

