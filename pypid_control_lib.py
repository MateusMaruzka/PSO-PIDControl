#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 20:03:57 2019

@author: Prof. Daniel Cavalcanti Jeronymo - UTFPR
@author: Mateus Maruzka Roncaglio - UTFPR
"""
import numpy as np
import math
import sys

import scipy.signal
import matplotlib.pyplot as plt

def simc(ts, tau):
    theta = ts
    t = tau
    k = 1
    kp = (2*t + theta)/(3*theta*k)
    ti = min(t+theta/2, 8*theta)
    
    return np.array([[kp, kp/ti]])

def isimc(ts,tau):
    theta = ts
    k=1
    t = tau# da func de transf
    tc = theta
    ti = min(t+theta/3, 4*(tc+theta))
    kc = (t+theta/3)/(k*(tc + theta))
    return np.array([[kc, kc/ti]])
    
def zn(ts,tau):
    Kp_zn = 0.9*tau/ts
    Ti_zn = 3*ts
    return np.array([[Kp_zn, Kp_zn/Ti_zn]])


def imc(L, T, k = 1):#
    
    theta = L
    tau = T

    l = 2*(tau+theta)/3
    #l = theta #skogestad 03
    kp = tau/(k*(l+theta))
    ti = min(tau,4*(l+theta)) 

    return np.array([[kp, kp/ti]])
    

def cc(ts,tau):
    theta = ts
    t = tau
    k=1
    kp = (t/theta)*(0.9 + theta/(12*t))/k
    ti = theta*(30+3/t)/(13+8/t)

    return np.array([kp, 1/ti])



def ise3(x):
    
    # Ainda pode ocorrer overflow
    B = np.any(x>np.sqrt(sys.float_info.max-1), axis = 1, keepdims = 1) 
    # print(B)
    e = np.sum(x**2, axis = 1, where = ~B)
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

def ise(e):
    
    # return ise2(e)
    return np.sum(e**2,axis=1)

def iae(e):
    return np.sum(np.abs(e),axis=1)

def itae(e):
    k = np.arange(len(e[0]))
    return np.dot(np.abs(e),k)

def tvc(u):
    return ise(u[...,1:-1] - u[...,0:-2])

def step_info(t,yout): 
    #t = iter(t1)
    os = (yout.max()/yout[-1]-1)*100
    tr = t[next(i for i in range(0,len(yout)-1) if yout[i]>yout[-1]*.90)]-t[0]
    #print ("OS: %f%s"%((yout.max()/yout[-1]-1)*100,'%'))
    #print ("Tr %f"%(t[next(i for i in range(0,len(yout)-1) if yout[i]>yout[-1]*.90)]-t[0]))
    A = abs(yout - 1) < 0.02 # ts

    tp = np.argmax(yout)
    # tp = tp*t[tp]
    try:
        ts = t[A][0]
    except IndexError:
        ts = 99
   # print("Ts %f"%t[A][0])
    return os, tr, ts, tp

def step(P, ts, tf, atraso = 0):
    
    Pd = P.to_discrete(ts)
    
    B = Pd.num                  # zeros
    A = Pd.den                  # poles
    nb = len(B) - 1             # number of zeros
    na = len(A) - 1             # number of poles
        
    slack = np.amax([na, nb]) + 1 + atraso # slack for negative time indexing of arrays
    kend = math.ceil(tf/ts) + 1   # end of simulation in discrete time
    kmax = kend + slack           # total simulation array size
    t = np.arange(0, tf + ts, ts)

    
    y = np.zeros(kmax)
    u = np.ones(kmax)
    u[0:atraso] = 0


    # Simulate
    for k in range(slack, kmax):
        y[k] = np.dot(B, u[k-1-atraso:k-1-atraso-(nb+1):-1]) - np.dot(A[1:], y[k-1:k-1-na:-1])
       
    y = y[slack:]
    t = np.arange(0, tf + ts, ts)

    return y, t


def picontrol(P, ts, tf, vetor_ganhos, num_controladores, atraso = 0):
    
    Pd = P.to_discrete(ts)
    
    B = Pd.num                  # zeros
    A = Pd.den                  # poles
    nb = len(B) - 1             # number of zeros
    na = len(A) - 1             # number of poles
        
    slack = np.amax([na, nb]) + 1 + atraso # slack for negative time indexing of arrays
    kend = math.ceil(tf/ts) + 1   # end of simulation in discrete time
    kmax = kend + slack           # total simulation array size
    t = np.arange(0, tf + ts, ts)

    
    y =  np.zeros([kmax, num_controladores])
    yD = np.zeros([kmax, num_controladores])
    u =  np.zeros([kmax, num_controladores])
    e =  np.zeros([kmax, num_controladores])
    r =  np.ones([kmax, num_controladores])

    yD[kmax//2:] = 0
 
    kp = vetor_ganhos[:,0]
    ki = vetor_ganhos[:,1]

    # Simulate
    for k in range(slack, kmax):
        y[k] = np.dot(B, u[k-1-atraso:k-1-atraso-(nb+1):-1]) - np.dot(A[1:], y[k-1:k-1-na:-1])
        
        y[k] = y[k] + yD[k]
       
        # error
        e[k] = r[k] - y[k]
        
        # PI control discretized by backwards differences
        du = kp*(e[k] - e[k-1])  + ki*e[k]*ts
        u[k] = u[k-1] + du 
      
        
  
    y = y[slack:]
    e = e[slack:]
    u = u[slack:]
    r = r[slack:]
        
    return y.T, e.T, u.T, r.T, t



def main():
    
    
    Ts = 0.1
    
    # P = scipy.signal.TransferFunction([1], [1, 4, 6, 4, 1])
    P = scipy.signal.TransferFunction([1], [2.3335711, 1])

    
    y,e,u,r,t = picontrol(P, Ts, 200, vetor_ganhos=np.array([[15.02626755, 19.57202986]]), num_controladores=1, atraso=18)
    

    print(np.sum(e**2, axis=1))
    plt.plot(t,e[0] ,'k.')
    # plt.plot(t,y[1] ,'k.')

    


if __name__ == "__main__":
    main()


