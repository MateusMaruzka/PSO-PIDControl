#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 10:27:39 2019

@author: Mateus Maruzka 
"""
import numpy as np
import scipy.signal
import math

from pypso_lib import atualizaVel, atualizaFitness, coefInercial, pso, limitaRegiao
from pypid_control_lib import ise


def picontrol(P, ts, tf,x, T_ENXAME):
    
    # Parametros da simulacao do controle pi
    #P = scipy.signal.TransferFunction([1], [10, 1])
    #P = scipy.signal.TransferFunction([1], [1, 9,23,15]) 
    Pd = P.to_discrete(ts)
    
    B = Pd.num                  # zeros
    A = Pd.den                  # poles
    nb = len(B) - 1             # number of zeros
    na = len(A) - 1             # number of poles
        
    slack = np.amax([na, nb]) + 1 # slack for negative time indexing of arrays
    kend = math.ceil(tf/ts) + 1   # end of simulation in discrete time
    kmax = kend + slack           # total simulation array size
    
    y = np.zeros([kmax, T_ENXAME])
    u = np.zeros([kmax, T_ENXAME])
    e = np.zeros([kmax, T_ENXAME])
    r = 1*np.ones([kmax, T_ENXAME])
    
    kp = x[:,0]
    ki = x[:,1]
    
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
        
    return e.T, u.T

def func_fitness(x,P,ts,tf,LAMBDA):

    mE, mDU = picontrol(P,ts,tf,x,len(x)) # testar
    mDU = mDU[...,1:-1] - mDU[...,0:-2]
    f = ise(mE) + LAMBDA*ise(mDU) # MULTIOBJETIVO (LQR)
    
    return f
    
#def func_pso(
#        T_ENXAME = 50,    
#        Imax = 20, #Numero maximo de iteraçoes para calculo do coef inercial (w)
#        Wmax = 0.8,#coef inercial max
#        Wmin = 0.4,
#        alfa = 1,
#        c1 = 2,
#        c2 = 2#coef inercial min   
#   ):
def func_pso(X):
    
    X = X[0]
    print("oi")
    print(X)
    DIM = len(X)     
    T_ENXAME = int(X[0])
    Imax = X[1]
    Wmax = X[2]
    Wmin = X[3]
    alfa = X[4]
    c1 = X[5]
    c2 = X[6]
    
    P = scipy.signal.TransferFunction([1], [1, 4, 6, 4, 1])

    ts = 0.1
    tf = 10
    
#------ Inicia váriaveis -------------    
    x = alfa*np.random.randn(T_ENXAME,DIM)

    v = np.random.randn(T_ENXAME,DIM)    
    pBest = x #pbest recebe o valor de x por ser a única posição conhecida da partic
    fitPbest = np.inf*np.ones(len(x)) 
    fitIter = []
    #k = 2/np.abs(2-(c1+c2)-np.sqrt((c1+c2)**2-4*(c1+c2)))
    #def atualizaFitness(func_fitness, posAtual, fitpBest, pbest, **args):
    # def func_fitness(x,P,ts,tf,LAMBDA):
    atualizaFitness(func_fitness,x, fitPbest, pBest, P=P,ts=ts,tf=tf,LAMBDA=0) # atualiza fitness atual e pBest 
    
    gb = np.argmin(fitPbest)
    fitIter.append(fitPbest[gb])
    for j in range(10):
        
        w = coefInercial(Wmin,Wmax,j,Imax)
        # def atualizaVel(x,v,pbest, gbest, num_particulas, w, c1 = 2, c2 = 2):
        v = atualizaVel(x,v,pBest,pBest[gb],T_ENXAME,w,c2,c1)  #wIter é um vetor com os valores de w a cada iteraçao
        # limitaVelocidade(v, 5000)
        x = x + v
       # limitaRegiao(x, 500)
        atualizaFitness(func_fitness,x, fitPbest, pBest, P=P,ts=ts,tf=tf,LAMBDA=0) # atualiza fitness atual e pBest 
        gb = np.argmin(fitPbest)
        fitIter.append(fitPbest[gb])
          
    #return pBest[gb], fitIter[]
    return fitIter.pop()


    #f = func_fitness(posAtual)
    #A = f < fitpBest # vetor de decisao
    #fitpBest[A] = f[A]
    #pbest[A] = posAtual[A]
    
def mypso(fObj,T_ENXAME, DIM):
    
    x = 10*np.random.rand(T_ENXAME,DIM)
    v = np.random.randn(T_ENXAME,DIM)
    Wmin = 0.1
    Wmax = 0.9
    c1 = 2
    c2 = 2
    pBest = x
    fitPbest = np.inf*np.ones(len(x)) 
    i = 0
    
    while i < 10:
        
        w = coefInercial(Wmin,Wmax,i,50)
        atualizaFitness(fObj, x,fitPbest,pBest) # atualiza fitness atual e pBest 
        gb = np.argmin(fitPbest) # gb = indice da particula com a melhor posiçao
        v = atualizaVel(x,v,pBest,pBest[gb], T_ENXAME,w,c2,c1)  #wIter é um vetor com os valores de w a cada iteraçao
        #  elitismo = x[gb]
        x = x + v
        limitaRegiao(x, 100)
        i = i + 1
        
    return pBest[gb]

def main():
    print(mypso(func_pso, 1, 7))   


if __name__ == "__main__":
    main()
