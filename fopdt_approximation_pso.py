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
import matplotlib.pyplot as plt



# Parametros da simulacao do controle pi
#P = scipy.signal.TransferFunction([1], [1, 9,23,15])
T_ENXAME = 50

def limitaRegiao(x, x_lim):
    A = x > x_lim
    x[A] = x_lim
    A = x < 0
    x[A] = 0

def ise(e):
    return np.sum(e**2,axis=1)

def iae(e):
    return np.sum(np.abs(e),axis=1)

def itae(e):
    k = np.arange(len(e))
    return np.dot(np.abs(e),k)


def compara_sinais(x,y1,t1):
    fit = np.zeros(len(x))
    for i in range(len(x)):
        L = x[i][0]
        T = x[i][1]
        G = scipy.signal.TransferFunction(np.polymul([1], [-L/2,1]), np.polymul([T,1], [L/2, 1]))

        t,y2 = scipy.signal.step(G,T =t1)
        fit[i] = np.sum((y2-y1)**2)
    return fit

def atualizaFitness(posAtual,fitpBest,pbest):
    
    P = scipy.signal.TransferFunction([1], [1, 4, 6, 4, 1])

    t1,y1 = scipy.signal.step(P)
    f = compara_sinais(posAtual,y1,t1)
    A = f < fitpBest # vetor de decisao
    fitpBest[A] = f[A]
    pbest[A] = posAtual[A]
    
    
    
def atualizaVel(x,v,pbest, gbest):
    
   # new_vel = np.zeros([T_ENXAME,DIM])
    c1 = 2
    r1 = np.random.rand() # entre 0 e 1
    c2 = 2
    r2 = np.random.rand() # entre 0 e 1
   
#    if j < Imax:
#        w = Wmax - (Wmax - Wmin)*j/Imax
#    else:
#        w = Wmin
#    www.append(w)

    #vel(t+1) =  w*vel(t) + c1*r1*(Pbest(t) - Pos(t))+ c2*r2*(gbest(t) - Pos(t))
    return 0.8*v + c1*r1*(pbest - x) + c2*r2*(np.tile(gbest,[T_ENXAME,1])-x)
   

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
        v = atualizaVel(x,v,pBest,pBest[gb])  #wIter é um vetor com os valores de w a cada iteraçao
        #  elitismo = x[gb]
        x = x + v
        limitaRegiao(x, 100)
    
    return pBest[gb]
    
def step_info(t,yout): 
    #t = iter(t1)
    
    print ("OS: %f%s"%((yout.max()/yout[-1]-1)*100,'%'))
    print ("Tr %f"%(t[next(i for i in range(0,len(yout)-1) if yout[i]>yout[-1]*.90)]-t[0]))
    A = abs(yout - 1) < 0.02 # ts
    print("Ts %f"%t[A][0])
    
   
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
    
    colors = ['b','g','r','m','k',]
    
    if len(x[0]) != len(y[0]):
        print("Algo de errado não está certo")
    else:   
        lines = []
        t = iter(colors)
        for i,j in zip(x, y):
            line, = plt.step(i, j, '-',  color=next(t), linewidth=1.2)
            lines.append(line)
    #ax1.plot(x, y1, color=[0.5, 0.5, 0.5], linewidth=1.2)
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    
    # plt.title(title, loc='center')  
    
    # melhor lugar de legenda
    #plt.legend()
    plt.legend(tuple(lines), tuple(legends), loc = 'best')
    # plt.legend((line6, line1, line2, line3, line4, line5, line6), (legend6, legend1, legend2, legend3, legend4, legend5), loc='best')
    plt.show()
    
  # plt.savefig("./stoi_result.svg")
    
   # plt.savefig("./stoi_result.jpg")
    


def main():
   gBest = pso(50, 2)
   print(gBest)
   L = gBest[0]
   T = gBest[1]
   t1,y1 = scipy.signal.step(scipy.signal.TransferFunction(np.polymul([1], [-L/2,1]), np.polymul([T,1], [L/2, 1])), N=1000)
   t2,y2 = scipy.signal.step(scipy.signal.TransferFunction([1], [1, 4, 6, 4, 1]), T = t1)

   step_info(t1,y1)
   step_info(t2,y2)
   
   maruzka_plot((t1,t2),(y1,y2), xlabel = "t (s)", ylabel = "Saída",legends=["FOPDT","4º Ordem"])
if __name__ == "__main__":
    main()
         

