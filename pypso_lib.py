#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 20:18:41 2019

@author: Mateus Maruzka Roncaglio
@author: Prof Daniel Cavalcanti Jeronymo 
@author: Wesley Klewerton 
"""
import numpy as np
import scipy.signal

def limitaRegiao(x, x_lim):
    A = x > x_lim
    x[A] = x_lim
    A = x < 0
    x[A] = 0


def atualizaVel(x,v,pbest, gbest, num_particulas, func_coef_inercial, c1 = 2, c2 = 2):
    
   # new_vel = np.zeros([T_ENXAME,DIM])
    r1 = np.random.rand() # entre 0 e 1
    r2 = np.random.rand() # entre 0 e 1

    w = func_coef_inercial(atualizaVel.iteracao)
    return w*v + c1*r1*(pbest - x) + c2*r2*(np.tile(gbest,[num_particulas,1])-x)


def atualizaFitness(func_fitness, posAtual, fitpBest, pbest):
    
    # Ts = 0.1
    # P = scipy.signal.TransferFunction([1], [1, 4, 6, 4, 1])
    # Pd = P.to_discrete(Ts)
    #y1,t1 = step(Pd.num, Pd.den, Ts, 20)
    t1,y1 = scipy.signal.step(P)
    # f = compara_sinais(posAtual,y1,t1, Ts)
    f = func_fitness(posAtual)
    A = f < fitpBest # vetor de decisao
    fitpBest[A] = f[A]
    pbest[A] = posAtual[A]
    



def pso(T_ENXAME, DIM):
    
    x = 11*np.random.rand(T_ENXAME,DIM)
    v = np.random.randn(T_ENXAME,DIM)
    
    #pbest recebe o valor de x por ser a única posição conhecida da partic
    pBest = x
    fitPbest = np.inf*np.ones(len(x)) 
    i = 0
    while i < 100:
        i = i + 1
        atualizaFitness(x,fitPbest,pBest) # atualiza fitness atual e pBest 
        gb = np.argmin(fitPbest) # gb = indice da particula com a melhor posiçao
        v = atualizaVel(x,v,pBest,pBest[gb], T_ENXAME)  #wIter é um vetor com os valores de w a cada iteraçao
        #  elitismo = x[gb]
        x = x + v
        limitaRegiao(x, 100)
    
    return pBest[gb]
