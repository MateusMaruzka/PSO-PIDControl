# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 17:47:34 2018

@author: Maruzketi
"""

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(0)

def schaffer_f6(x):
    
    num = np.sin(np.degrees(np.sqrt(np.sum(x**2,axis=1))))**2 - 0.5
    den = (1 + 0.001*(np.sum(x**2,axis=1)))
    
    return (0.5+ num/den)    
    

def sphereFunction_v2(x):
    return np.sum(x**2,axis=1)



def imprimeDados(fitIter,gBest):
    print('Melhor solucao:')
    print(gBest)
    print('Fitness da melhor solucao')
    print(fitIter[-1])
    plt.figure()
    plt.plot(fitIter)
   # plt.figure()
   # plt.plot((wIter))

def inicia_variaveis(DIM, T_ENXAME):
    sigma1 = 5
    x = sigma1*np.random.randn(T_ENXAME,DIM)
    v = np.random.randn(T_ENXAME,DIM)     
    pBest = x
    fit_pBest = np.inf*np.ones(len(x))
    return x,v,pBest, fit_pBest

def atualizaFitness(posAtual,fitpBest,pbest):
    
    y = sphereFunction_v2(posAtual)
    #y = sphereFunction_v2(posAtual)
    A = y < fitpBest # vetor de decisao
    fitpBest[A] = y[A]
    pbest[A] = posAtual[A]



def atualizaVel(x, v, pbest, gbest, j, T_ENXAME, DIM, Wmax=0.9, Wmin=0.05, Imax = 100):
    new_vel = np.zeros([T_ENXAME,DIM])
    c1 = 1
    r1 = np.random.rand(T_ENXAME,1) 
    c2 = 1
    r2 = np.random.rand(T_ENXAME,1) 
          
    if j < Imax:
        w = Wmax - (Wmax - Wmin)*j/Imax
    else:
        w = Wmin
        
   # sl_factor = np.sum(pbest,axis=0)/len(pbest) # share-learning factor  #+ c3*r3*(np.tile(sl_factor,[T_ENXAME,1])-x) 
    new_vel = w*v + c1*r1*(pbest - x) + c2*r2*(np.tile(gbest,[T_ENXAME,1])-x) 
   
    
    return new_vel



def pso(Iter_max, DIM,T_ENXAME):
    """
    PSO padrão:
        Retorna a melhor fitness a cada iteração através de fitIter
        Retorna gBest
    """
    x,v,pBest,fitPbest = inicia_variaveis(DIM, T_ENXAME)
    fitIter = []
    j = 0
    while j < Iter_max:
        j+=1    
        atualizaFitness(x,fitPbest,pBest) # atualiza fitness atual e pBest 
        gb = np.argmin(fitPbest) # gb = indice da particula com a melhor posiçao
        v = atualizaVel(x,v,pBest,pBest[gb],j,T_ENXAME,DIM)  #wIter é um vetor com os valores de w a cada iteraçao
        x = x + v
        fitIter.append(fitPbest[gb])
    
    return fitIter, pBest[gb]



def main():
    DIM = 5  #dimensoes do problema
    T_ENXAME = 30  #tamanho do enxame
    Iter_max = 100
    gBest = 0;
    fitIter = []
    
    fitIter, gBest = bbpso(Iter_max, DIM, T_ENXAME)
    imprimeDados(fitIter, gBest)
    fitIter, gBest = pso(Iter_max, DIM, T_ENXAME)
    imprimeDados(fitIter, gBest)


if __name__ == "__main__":
    main()


