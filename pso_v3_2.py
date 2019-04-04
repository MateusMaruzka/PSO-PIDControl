# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 16:11:04 2017

@author: Maruzka
"""

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(2)

def schaffer_f6(x):
    
    num = np.sin(np.degrees(np.sqrt(np.sum(x**2,axis=1))))**2 - 0.5
    den = (1 + 0.001*(np.sum(x**2,axis=1)))
    
    return (0.5+ num/den)    
    

def sphereFunction_v2(x):
    return np.sum(x**2,axis=1)


def limitaRegiao(x, x_lim):
    A = x > x_lim
    x[A] = x_lim

def limitaVelocidade(v, v_lim):
    A = v > v_lim
    v[A] = v_lim

def atualizaFitness(posAtual,fitpBest,pbest, function = 'sphere'):
    
    #y = schaffer_f6(posAtual)
    if function == 'sphere':
        y = sphereFunction_v2(posAtual)
    elif function == 'schaffer':
        y = schaffer_f6(posAtual)
        
    A = y < fitpBest # vetor de decisao
    fitpBest[A] = y[A]
    pbest[A] = posAtual[A]


def atualizaVel(x,v,pbest,gbest,w,j,T_ENXAME,DIM,Wmax=0.9,Wmin=0.1,Imax = 50):
    new_vel = np.zeros([T_ENXAME,DIM])
    c1 = 1
    r1 = np.random.rand(T_ENXAME,1) 
    c2 = 1
    r2 = np.random.rand(T_ENXAME,1) 
          
    if j < Imax:
        w = Wmax - (Wmax - Wmin)*j/Imax
    else:
        w = Wmin
        
   # sl_factor = np.sum(pbest,axis=0)/len(pbest) # share-learning factor
    new_vel = w*v + c1*r1*(pbest - x) + c2*r2*(np.tile(gbest,[T_ENXAME,1])-x) 
    #+ c3*r3*(np.tile(sl_factor,[T_ENXAME,1])-x) 
    
    return new_vel

def imprimeDados(wIter, fitIter, fitPbest, pBest):
    gb = np.argmin(fitPbest)
    print('Melhor solucao:')
    print(pBest[gb])
    print('Fitness da melhor solucao')
    print(fitPbest[gb])
    plt.figure()
    plt.plot(fitIter)
    #plt.figure()
    #plt.plot((wIter))

    

def main():
    DIM = 3 #dimensoes do problema
    T_ENXAME = 50 #tamanho do enxame
    sigma1 = 5
    
    fig, ax = plt.subplots(2,1)
  
    for i in range(3):
        x = sigma1*np.random.randn(T_ENXAME,DIM)
        v = np.random.randn(T_ENXAME,DIM)     
        pBest = x
        fitPbest = np.inf*np.ones(len(x))
        fitIter = []
        wIter = []
        w = [0.9, 0.5, 0.2]
        j = 0
        while j < 100:
            j+=1    
            atualizaFitness(x,fitPbest,pBest) # atualiza fitness atual e pBest 
            gb = np.argmin(fitPbest) # gb = indice da particula com a melhor 
            fitIter.append(fitPbest[gb])
            v = atualizaVel(x,v,pBest,pBest[gb],wIter,j,T_ENXAME,DIM, Wmax = w[i])  #wIter é um vetor com os valores de w a cada iteraçao
            #limitaVelocidade(v,1000)
            x = x + v
            limitaRegiao(x,100)
            
        ax[0].plot(np.log10(fitIter), label = 'W=')
        #imprimeDados([], fitIter, fitPbest, pBest)
    plt.ylabel('Fitness')
    plt.xlabel('Iterações')
    plt.title('Convergência Função Esfera')
    ax[0].legend(w)
    
    
    for i in range(3):
        x = sigma1*np.random.randn(T_ENXAME,DIM)
        v = np.random.randn(T_ENXAME,DIM)     
        pBest = x
        fitPbest = np.inf*np.ones(len(x))
        fitIter = []
        wIter = []
        w = [0.9, 0.5, 0.1]
        j = 0
        while j < 100:
            j+=1    
            atualizaFitness(x,fitPbest,pBest, function='schaffer') # atualiza fitness atual e pBest 
            gb = np.argmin(fitPbest) # gb = indice da particula com a melhor 
            fitIter.append(fitPbest[gb])
            v = atualizaVel(x,v,pBest,pBest[gb],wIter,j,T_ENXAME,DIM, Wmax = w[i])  #wIter é um vetor com os valores de w a cada iteraçao
            #limitaVelocidade(v,1000)
            x = x + v
            limitaRegiao(x,100)
            
        ax[1].plot(np.log10(fitIter), label = 'W=')

    plt.ylabel('Fitness')
    plt.xlabel('Iterações')
    plt.title('Convergência Função Esfera')
    ax[1].legend(w)
    

if __name__ == "__main__":
    main()


