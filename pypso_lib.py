#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 20:18:41 2019

@author: Mateus Maruzka Roncaglio
@author: Prof Daniel Cavalcanti Jeronymo 
@author: Wesley Klewerton 
"""
import numpy as np
#import scipy.signal

def coefInercial(Wmin, Wmax, i, imax):
    return (Wmax - (Wmax - Wmin)*i/imax)


def limitaRegiao(x, x_lim):
    A = x > x_lim
    x[A] = x_lim
    A = x < 0
    x[A] = 0


def atualizaVel(x,v,pbest, gbest, num_particulas, w, c1 = 2, c2 = 2):
    
   # new_vel = np.zeros([T_ENXAME,DIM])
    r1 = np.random.rand() # entre 0 e 1
    r2 = np.random.rand() # entre 0 e 1

   # w = func_coef_inercial(atualizaVel.iteracao)
    return w*v + c1*r1*(pbest - x) + c2*r2*(np.tile(gbest,[num_particulas,1])-x)


def atualizaFitness(func_fitness, posAtual, fitpBest, pbest, **args):
    
   # ff = np.vectorize(func_fitness)
   
    f = func_fitness(posAtual, **args)
    A = f < fitpBest # vetor de decisao
    fitpBest[A] = f[A]

    pbest[A] = posAtual[A]
    


def pso(fObj,T_ENXAME, DIM, iterMax = 50, _alfa = 30, _Wmin = 0.7,_Wmax = 0.9, _c1 = 2, _c2 = 2, **args):
    
    x = _alfa*np.random.rand(T_ENXAME,DIM)
    v = np.random.randn(T_ENXAME,DIM)
    pBest = x
    fitPbest = np.inf*np.ones(len(x)) 
    
    fitIter = []

    i = 0
    while i < iterMax:
        w = coefInercial(_Wmin,_Wmax,i,iterMax)
        atualizaFitness(fObj, x,fitPbest,pBest, **args) # atualiza fitness atual e pBest 
        gb = np.argmin(fitPbest) # gb = indice da particula com a melhor posiçao
        fitIter.append(fitPbest[gb])

        v = atualizaVel(x,v,pBest,pBest[gb], T_ENXAME,w,_c2,_c1)  #wIter é um vetor com os valores de w a cada iteraçao
        x = x + v
        
        limitaRegiao(x, 100)
        i = i + 1
        
    return pBest[gb],fitIter


def x_bbpso(pBest, fit_pBest, DIM, T_ENXAME):
    
    #np.tile(pBest[np.argmin(fit_pBest)],[T_ENXAME,1]) # Cria uma matriz com GBEST em cada linha
    mean = 0.5*(pBest + np.tile(pBest[np.argmin(fit_pBest)],[T_ENXAME,1]))
    std = np.abs(pBest - np.tile(pBest[np.argmin(fit_pBest)],[T_ENXAME,1]))
    
    return std*np.random.randn(T_ENXAME,DIM) + mean

def bbpso(P,ts,tf,T_ENXAME = 50):
    """
    Bare bone PSO
    """
    
    DIM = 2 #dimensoes do problema    
#------ Inicia váriaveis -------------    
    x = 2*np.random.randn(T_ENXAME,DIM)
    pBest = x #pbest recebe o valor de x por ser a única posição conhecida da partic
    fitPbest = np.inf*np.ones(len(x)) 
    fitIter = []
    atualizaFitness(P,ts,tf,fitPbest,pBest,x,0) # atualiza fitness atual e pBest 
    gb = np.argmin(fitPbest)
    fitIter.append(fitPbest[gb])
    for j in range(30):
        
        x = x_bbpso(pBest, fitPbest, DIM, T_ENXAME)
       # limitaRegiao(x, 5000)
        atualizaFitness(P,ts,tf,fitPbest,pBest,x,0) # atualiza fitness atual e pBest 
        gb = np.argmin(fitPbest)
        fitIter.append(fitPbest[gb])

    return pBest[gb], fitIter
    #return fitIter