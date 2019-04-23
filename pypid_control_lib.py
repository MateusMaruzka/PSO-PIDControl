#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 20:03:57 2019

@author: Prof. Daniel Cavalcanti Jeronymo - UTFPR
@author: Mateus Maruzka Roncaglio - UTFPR
"""
import numpy as np
import math

def imc(ts, tau):#
    
    theta = ts
    k=1
    l = 2*(tau+theta)/3
    l = theta #skogestad 03
    kp = tau/(k*(l+theta))
    ti = min(tau,4*(l+theta)) 

    return([kp, 1/ti])
    
def simc(ts, tau):
    theta = ts
    t = tau
    k = 1
    kp = (2*t + theta)/(3*theta*k)
    ti = min(t+theta/2, 8*theta)
    
    return np.array([kp, 1/ti])

def isimc(ts,tau):
    theta = ts
    k=1
    t = tau# da func de transf
    tc = theta
    ti = min(t+theta/3, 4*(tc+theta))
    kc = (t+theta/3)/(k*(tc + theta))
    return np.array([kc, kc/ti])
    
def zn(ts,tau):
    Kp_zn = 0.9*tau/ts
    Ti_zn = 3*ts
    return np.array([Kp_zn, 1/Ti_zn])



def cc(ts,tau):
    theta = ts
    t = tau
    k=1
    kp = (t/theta)*(0.9 + theta/(12*t))/k
    ti = theta*(30+3/t)/(13+8/t)

    return np.array([kp, 1/ti])


def ise(e):
    return np.sum(e**2,axis=1)

def iae(e):
    return np.sum(np.abs(e),axis=1)

def itae(e):
    k = np.arange(len(e))
    return np.dot(np.abs(e),k)

def step_info(t,yout): 
    #t = iter(t1)
    
    print ("OS: %f%s"%((yout.max()/yout[-1]-1)*100,'%'))
    print ("Tr %f"%(t[next(i for i in range(0,len(yout)-1) if yout[i]>yout[-1]*.90)]-t[0]))
    A = abs(yout - 1) < 0.02 # ts
    print("Ts %f"%t[A][0])
    

def step(b, a, Ts = 0.1, tf = 1):
    
    # Pd = P.to_discrete(Ts)  
    # B = Pd.num                  # zeros
    # A = Pd.den                  # poles
    B = b
    A = a
    nb = len(B) - 1             # number of zeros
    na = len(A) - 1             # number of poles

    # Simulation parameters
    # tf = 0.05
    
    slack = np.amax([na, nb]) + 1 # slack for negative time indexing of arrays
    kend = math.ceil(tf/Ts) + 1   # end of simulation in discrete time
    kmax = kend + slack           # total simulation array size
    
    y = np.zeros(kmax)
    u = np.ones(kmax)
    # u[0:slack] = 0
    
    # Simulate
    for k in range(slack, kmax-1):
        #print(len(A[1:]))
        #print(len(y[k-1:k-1-na:-1]))
        y[k] = np.dot(B, np.ones(len(b))) - np.dot(A[1:], y[k-1:k-1-na:-1])

        #y[k] = np.dot(B, u[k-1:k-1-(nb+1):-1]) - np.dot(A[1:], y[k-1:k-1-na:-1])
        
    # ignore empty slack
    y = y[slack:]
    u = u[slack:]
    t = np.arange(0, tf + Ts, Ts)
    
    return y,t

def picontrol(P, ts, tf, vetor_ganhos, num_controladores):
    
    Pd = P.to_discrete(ts)
    
    B = Pd.num                  # zeros
    A = Pd.den                  # poles
    nb = len(B) - 1             # number of zeros
    na = len(A) - 1             # number of poles
        
    slack = np.amax([na, nb]) + 1 # slack for negative time indexing of arrays
    kend = math.ceil(tf/ts) + 1   # end of simulation in discrete time
    kmax = kend + slack           # total simulation array size
    
    y = np.zeros([kmax, num_controladores])
    u = np.zeros([kmax, num_controladores])
    e = np.zeros([kmax, num_controladores])
    r = 1*np.ones([kmax, num_controladores])
    
    kp = vetor_ganhos[:,0]
    ki = vetor_ganhos[:,1]

    # Simulate
    for k in range(slack, kmax):
        y[k] = np.dot(B, u[k-1:k-1-(nb+1):-1]) - np.dot(A[1:], y[k-1:k-1-na:-1])
        
        # error
        e[k] = r[k] - y[k]
        
        # PI control discretized by backwards differences
        du = kp*(e[k] - e[k-1])  + ki*e[k]*ts
        u[k] = u[k-1] + du 
         
        # SATURACAO
        # u[k] = min(max(0, u[k]), 2)
        # gg = u[k] > 10
        # u[k,gg] = 10
        
    print(slack)
    print("resp")
    print(y.T)
    #print(e.T)
    #print(u.T)
    return y.T[...,slack:], e.T[...,slack:], u.T[...,slack:]


# Plot time response
"""

h(z) = y(z)/u(z) = pol(a)/pol(b)

y(z)*pol(b) = pol(a)*u(z)

a(n)*z^(n) + a(n-1)*z^(n-1) +a(n-2)*z^(n-2) ...

b(n)*z^(n) + b(n-1)*z^(n-1) +b(n-2)*z^(n-2) ...

pol(a)*u(z) = [z^(m) + a(1)*z^(m-1) +a(2)*z^(m-2)  + An ...]*[U(z)]

pol(b)*y(z) = [z^(n) + b(1)*z^(n-1) +b(2)*z^(n-2) + ... +Bn ]*[Y(z)]


[z^(n) + b(1)*z^(n-1) +b(2)*z^(n-2) + ... +Bn ]*[Y(z)] = [z^(m) + a(1)*z^(m-1) +a(2)*z^(m-2)  + An ...]*[U(z)]

y[n] + b1*y[n-1] +b2*y[n-2] + ... + bn = u[z] + a1*y[m-1] + a2*z[m-2] + ... + an

y[n] = u[z] + a1*y[m-1] + a2*z[m-2] + ... + an - b1*y[n-1] +b2*y[n-2] + ... + bn

"""

