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
import types


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
    
def zn(atraso,tau, ctrl):
    
    if ctrl == "PI":
        
        Kp_zn = 0.9*tau/atraso
        Ti_zn = atraso/0.3
        
        return np.array([[Kp_zn, Kp_zn/Ti_zn]])
    
    elif ctrl == "PID":
        
        Kp_zn = 1.2*tau/atraso
        Ti_zn = 2*atraso
        Kd_zn = 0.5*atraso
        
        return np.array([[Kp_zn, Kp_zn/Ti_zn, Kd_zn]])
    
    else:
        raise AttributeError("Atributo ctrl deve ser 'PI' ou 'PID' ")
    



def imc(atraso, tau, tauc, K = 1):
    
    kc = (1/K)*(tau/(tauc+atraso))
    ti = tau
    
    return np.array([[kc, kc/ti]])


def imc_old(L, T, k = 1):#
    
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

    return np.sum((u[...,1:-1] - u[...,0:-2])**2)

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
    return os, tr, ts

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


def picontrol(P, ts, tf, vetor_ganhos, num_controladores, atraso = 0, d = 0, pid_discrete_function=0):
    

    if isinstance(pid_discrete_function, types.FunctionType):
        picontrol.func = pid_discrete_function
    else:
        
        picontrol.func = lambda e,u,k, kp, ki: kp*(e[k] - e[k-1])  + ki*e[k]*ts + u[k-1]
        # def picontrol_function(e, u, k, kp, ki):
        #     return kp*(e[k] - e[k-1])  + ki*e[k]*ts + u[k-1]
        
        # picontrol.func = picontrol_function
    
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

    yD[kmax//2:] = d
 
    kp = vetor_ganhos[:,0]
    ki = vetor_ganhos[:,1]

    # print(kp, ki)
    # Simulate
    for k in range(slack, kmax):
        y[k] = np.dot(B, u[k-1-atraso:k-1-atraso-(nb+1):-1]) - np.dot(A[1:], y[k-1:k-1-na:-1])
        
        y[k] = y[k] + yD[k]
       
        # error
        e[k] = r[k] - y[k]
        
        
        u[k] = picontrol.func(e, u,k, kp, ki)
        # PI control discretized by backwards differences
        # u[k] - u[k-1] = (kp+ki*ts)e[k] - kp*e[k-1]
        #  z - 1 = (kp+ki*ts)z - kp
        # Cpid = [(kp+ki*ts)z - kp]/(z - 1)
        
        
        # u[k] = picontrol.discrete_pid_function(e, u)
        # du = kp*(e[k] - e[k-1])  + ki*e[k]*ts 
        # u[k] = u[k-1] + du 
        
  
    y = y[slack:]
    e = e[slack:]
    u = u[slack:]
    r = r[slack:]
        
    return y.T, e.T, u.T, r.T, t



def main():
    
    
    Ts = 0.1
    
    P = scipy.signal.TransferFunction([1], [1, 4, 6, 4, 1])
    # P = scipy.signal.TransferFunction([1], [2.3335711, 1])

    
    y,e,u,r,t = picontrol(P, Ts, 100, vetor_ganhos=np.array([[0.5, .1]]),\
                          num_controladores=1, atraso=0, pid_discrete_function= lambda e,u,k,kp,ki: 1/k)
    

    print(np.sum(e**2, axis=1))
    plt.plot(t,y[0] ,'k.')
    # plt.plot(t,y[1] ,'k.')

    

        
        
if __name__ == "__main__":
    
    if isinstance(main, types.FunctionType):
        print("sucesso")
    main()
    


