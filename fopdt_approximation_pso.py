# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 18:02:06 2018
@author: Maruzka
"""

"""
Otimização por enxame de partículas para encontrar L e T de um FOPDT

F(s) = exp(-Ls)/(Ts+1) 

"""


import numpy as np
import scipy.signal
import math
import matplotlib.pyplot as plt

from pypid_control_lib import step, step_info
from pypso_lib import limitaRegiao
# Parametros da simulacao do controle pi
#P = scipy.signal.TransferFunction([1], [1, 9,23,15])
#T_ENXAME = 50

def linearW(i,Imax, Wmax = 0.9, Wmin = 0.2):
   return Wmax - (Wmax - Wmin)*i/Imax
    
def compara_sinais(x,y1,t1, Ts):
    
    fit = np.zeros(len(x))
    for i in range(len(x)):
        L = int(x[i][0])
        x[i][0] = L
        T = x[i][1]
        # G = scipy.signal.TransferFunction(np.polymul([1], [-L/2,1]), np.polymul([T,1], [L/2, 1]))
        if T < 0.01:
            x[i][1] = 0.01
            T = 0.01
            
        G = scipy.signal.TransferFunction([1], [T, 1])
        # y2,t = step(Gd.num, Gd.den, Ts, 20)
        y2,t = step(G, Ts, tf = 10, atraso=L)
        # t,y2 = scipy.signal.step(G,T =t1)
        #fit[i] = np.dot((y2-y1)**2, np.arange(len(y2)))
        fit[i] = np.sum((y2-y1)**2)
        #t = np.arange(len(y1))
        #fit[i] = np.sum(np.dot(np.abs(y2-y1), t))
    return fit

def atualizaFitness(posAtual,fitpBest,pbest):
    
    Ts = 0.1

    P = scipy.signal.TransferFunction([1], [1, 4, 6, 4, 1])
    y1,t1 = step(P, Ts, tf = 10)

    f = compara_sinais(posAtual,y1,t1, Ts)
    A = f < fitpBest # vetor de decisao
    fitpBest[A] = f[A]
    pbest[A] = posAtual[A]
    
    
    
def atualizaVel(x,v,pbest, gbest, num_particulas, func_coef_inercial, c1 = 2.05, c2 = 2.05):
    atualizaVel.iteracoes = atualizaVel.iteracoes + 1

    r1 = np.random.rand() # entre 0 e 1
    r2 = np.random.rand() # entre 0 e 1
   
    w = func_coef_inercial(atualizaVel.iteracoes,100)
    
    return w*v*0.86 + 0.86*c1*r1*(pbest - x) + 0.86*c2*r2*(np.tile(gbest,[num_particulas,1])-x)
   

def pso(T_ENXAME, DIM):
    
    x = 11*np.random.rand(T_ENXAME,DIM)
    v = np.random.randn(T_ENXAME,DIM)
    
    #pbest recebe o valor de x por ser a única posição conhecida da partic
    pBest = x
    fitPbest = np.inf*np.ones(len(x)) 
    atualizaVel.iteracoes = 0
    i = 0
    while i < 100:
        i = i + 1
        atualizaFitness(x,fitPbest,pBest) # atualiza fitness atual e pBest 
        gb = np.argmin(fitPbest) # gb = indice da particula com a melhor posiçao
        v = atualizaVel(x,v,pBest,pBest[gb], T_ENXAME, linearW)  #wIter é um vetor com os valores de w a cada iteraçao
        x = x + v
        limitaRegiao(x, 100, 0)
    
    return pBest[gb]
   
def maruzka_plot(x,y, font_size = 12, figsize = (9,6), fontfamily = "Times New Roman", 
                 xlabel = "x label", ylabel = "y label", legends = ["legends1", "legends2"], 
                 title = "MyTitle"):
    
    #    figure generation      
    plt.rcParams['figure.figsize'] = figsize
    
    # tipo de fonte times new roman    
    plt.rcParams["font.family"] = fontfamily
    
    # tamanho da fonte da legenda
    plt.rcParams['legend.fontsize'] = font_size
    
    # cores
    colors = ['k','r','b','m','y',]
    
    if len(x[0]) != len(y[0]):
        print("Algo de errado não está certo")
    else:   
        lines = []
        t = iter(colors)
        for i,j in zip(x, y):
            line, = plt.plot(i, j, '-',  color=next(t), linewidth=1.2)
            lines.append(line)
    #ax1.plot(x, y1, color=[0.5, 0.5, 0.5], linewidth=1.2)
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    
    # plt.title(title, loc='center')  
    
    # melhor lugar de legenda
    plt.legend(tuple(lines), tuple(legends), loc = 'best')
    # plt.legend((line6, line1, line2, line3, line4, line5, line6), (legend6, legend1, legend2, legend3, legend4, legend5), loc='best')
    
    plt.savefig("./imagens/"+title+".svg")
    plt.show()


def main():
    
   gBest = pso(100, 2)
   print("gbest")
   print(gBest)
   L = gBest[0]
   T = gBest[1]
   print("")

   Ts = 0.1
   tf = 20
   # t = np.arange(0, tf + Ts, Ts)
   P1 = scipy.signal.TransferFunction([1], [T,1])

   y1,u = step(P1 ,Ts, 20, atraso=int(L))
   
   P2 = scipy.signal.TransferFunction([1], [1, 4, 6, 4, 1])
   
   y2,u = step(P2, Ts, 20)

   plt.plot(y2)
   plt.plot(y1)


if __name__ == "__main__":
    main()
         

