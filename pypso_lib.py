#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 20:18:41 2019

@author: Mateus Maruzka Roncaglio
@author: Prof Daniel Cavalcanti Jeronymo 
@author: Wesley Klewerton 
"""
import numpy as np
import matplotlib.pyplot as plt
#import scipy.signal

class constrictorFactorError(Exception):
        
        def __init__(self, *args):
            
            self.message = args if args else None
        
        

def coefInercial(Wmin, Wmax, i, imax):
    return (Wmax - (Wmax - Wmin)*i/imax)


def limitaRegiao(x, x_sup, x_inf):
    A = x > x_sup
    x[A] = x_sup
    A = x < x_inf
    x[A] = x_inf


def atualizaVel(x,v,pbest, gbest, num_particulas, w, c1 = 2, c2 = 2):
    
   # new_vel = np.zeros([T_ENXAME,DIM])
    r1 = np.random.random((v.shape)) # entre 0 e 1
    r2 = np.random.random((v.shape)) # entre 0 e 1

   # w = func_coef_inercial(atualizaVel.iteracao)
    return w*v + c1*r1*(pbest - x) + c2*r2*(gbest-x)


def atualizaFitness(func_fitness, posAtual, fitpBest, pbest, **args):
    
   # ff = np.vectorize(func_fitness)
   
    f = func_fitness(posAtual, **args)
    A = f < fitpBest # vetor de decisao
    fitpBest[A] = f[A]

    pbest[A] = posAtual[A]
    


def pso(fObj,T_ENXAME, DIM, iterMax = 100, var = 5, _Wmin = 0.1,_Wmax = 0.9, _c1 = 2.05, _c2 = 2.05, bondaries = None, **kwargs):
    
    
    constrictor_factor = lambda a,b : 2 / (np.abs(2-(a+b) - np.sqrt((a+b)**2 -4*(a+b))))

    try:
        cf = constrictor_factor(_c1, _c2)
        
        if(cf >= 1):
            raise constrictorFactorError("A soma de c1 e c2 deve ser maior que 4", (_c1,_c2))
            
    except constrictorFactorError as e:
        
        msg, _cf = e.args
        print('\033[1m'+'\033[33m'+ "Warning: "+ msg+'\033[1m', \
              "Fornecido: " + str(_cf), \
              "Corrigindo C1 e C2 para 2.05" + '\033[0;0m'\
              , sep="\n")
    
        cf = constrictor_factor(2.05, 2.05)
        _c1 = 2.05
        _c2 = 2.05
        
    
    
    fitIter = []

    
    x = var*np.random.rand(T_ENXAME,DIM)
    v = np.random.randn(T_ENXAME,DIM)
    pBest = np.copy(x)
    
    fitPbest = np.inf*np.ones(len(x)) 
    atualizaFitness(fObj, x,fitPbest,pBest, **kwargs)
    gb = np.argmin(fitPbest) 
    fitIter.append(fitPbest[gb])

    i = 0
    while i < iterMax:

        v = atualizaVel(x,v,pBest,pBest[gb], T_ENXAME,coefInercial(_Wmin,_Wmax,i,iterMax),_c2,_c1)  #wIter é um vetor com os valores de w a cada iteraçao
        x = x + cf*v
        limitaRegiao(x, x_sup = 100, x_inf = 0.01)
        
        atualizaFitness(fObj, x,fitPbest,pBest, **kwargs) # atualiza fitness atual e pBest 
        gb = np.argmin(fitPbest) # gb = indice da particula com a melhor posiçao
        fitIter.append(fitPbest[gb])
        
        i = i + 1
        
    return pBest[gb],fitIter


def esfera(x):
    return np.sum((x)**2, axis=1)

def main():
    
    # np.random.seed(0)
    
    gb, fit = pso(esfera,10, 2, 50)

    print(gb)
    
    plt.semilogx(fit)
    
if __name__ == "__main__":
    main()
    
    
    
"""


def x_bbpso(pBest, fit_pBest, DIM, T_ENXAME):
    
    #np.tile(pBest[np.argmin(fit_pBest)],[T_ENXAME,1]) # Cria uma matriz com GBEST em cada linha
    mean = 0.5*(pBest + np.tile(pBest[np.argmin(fit_pBest)],[T_ENXAME,1]))
    std = np.abs(pBest - np.tile(pBest[np.argmin(fit_pBest)],[T_ENXAME,1]))
    
    return std*np.random.randn(T_ENXAME,DIM) + mean

def bbpso(P,ts,tf,T_ENXAME = 50):
    
    #Bare bone PSO
    
    
    DIM = 2 #dimensoes do problema    
#------ Inicia váriaveis -------------    
    x = 2*np.random.randn(T_ENXAME,DIM)
    pBest = x #pbest recebe o valor de x por ser a única posição conhecida da partic
    fitPbest = np.inf*np.ones(len(x)) 
    fitIter = []
    atualizaFitness(P,ts,tf,fitPbest,pBest,x,0) # atualiza fitness atual e pBest 
    gb = np.argmin(fitPbest)
    fitIter.append(fitPbest[gb])
    for j in range(50):
        
        x = x_bbpso(pBest, fitPbest, DIM, T_ENXAME)
       # limitaRegiao(x, 5000)
        atualizaFitness(P,ts,tf,fitPbest,pBest,x,0) # atualiza fitness atual e pBest 
        gb = np.argmin(fitPbest)
        fitIter.append(fitPbest[gb])

    return pBest[gb], fitIter
    #return fitIter


"""