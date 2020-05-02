# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 15:17:11 2020

@author: maruzka
"""

import numpy as np
import sys

np.seterr(all='warn')
np.seterr(over='raise', under = 'raise')



def mySum(x):
    
    # Ainda pode ocorrer overflow
    B = np.any(x>(sys.float_info.max//2), axis = 1, keepdims = 1) 
    # print(B)
    e = np.sum(x, axis = 1, where = ~B)
    # print(e)
    e[np.ravel(B)] = np.inf
    
    return e

def ise2(e):
    
    ise = np.zeros(len(e))
    
    for idx, i in enumerate(e):
        
        try:
            ise[idx] = np.sum(i**2)
        except FloatingPointError:
            print("Overflow")
            ise[idx] = np.inf
            
    return ise

x = np.arange(5000, dtype=np.float32).reshape(2,2500)
print(mySum(x))
print(np.sum(x, axis = 1))
        


# e = np.arange(9000000)*1E300
# np.append(e,e)

# def mySum():
    
#     pass

# try:
#     a = np.sum(e**2, axis = 1)
#     print("int16", a)
#     print("real", 32000*3)    

# except FloatingPointError as e:

#     try:
#         print("try")
#         e = np.array(e, dtype=np.float128)
#         a = np.sum(e**2)
    
#     except AttributeError:
#         print("vish")
#         e = np.inf