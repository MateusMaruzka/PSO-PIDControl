"""
Mateus Maruzka

Otimização por enxame de partículas otimizando um controlador PID

- Coeficiente inercial com decaimento por uma função linear
- Velocidade limite das partículas
- Restrição das partículas no espaço de busca
- Fator de constrição e valores de c1 e c2 calculados pela formula em shi
- Barebone pso

"""
from cycler import cycler
import scipy.signal
import numpy as np
import math
import matplotlib.pyplot as plt

from pypid_control_lib import imc, simc, isimc, zn, cc, iae, ise, itae, step_info
from pypso_lib import limitaRegiao

"""
padé approx

exp(-theta*s) = ((-theta/2)*s+1)/((theta/2)*s+1)
"""

def picontrol2(P,ts,tf,x):
    
    
    Pd = P.to_discrete(ts)    
    
    B = Pd.num                  # zeros
    A = Pd.den                  # poles
    nb = len(B) - 1             # number of zeros
    na = len(A) - 1             # number of poles
    
    slack = np.amax([na, nb]) + 1 # slack for negative time indexing of arrays
    kend = math.ceil(tf/ts) + 1   # end of simulation in discrete time
    kmax = kend + slack           # total simulation array size
    
    y = np.zeros(kmax)
    u = np.zeros(kmax)
    e = np.zeros(kmax)
    r = 1*np.ones(kmax)

    # Simulate
    kp = x[0]
    ki = x[1]
    
    for k in range(slack, kmax):
        y[k] = np.dot(B, u[k-1:k-1-(nb+1):-1]) - np.dot(A[1:], y[k-1:k-1-na:-1])
        
        # error
        e[k] = r[k] - y[k]
        
        # PI control discretized by backwards differences
        du = kp*(e[k] - e[k-1])  + ki*e[k]*ts
        u[k] = u[k-1] + du
         
        # SATURACAO
        #
        #u[k] = min(max(0, u[k]), 2)
        
#        gg = u[k] > 10
#        u[k,gg] = 10
 

    return y[slack:],u[slack:],e[slack:]
        
        
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
        #u[k] = min(max(0, u[k]), 2)
#        gg = u[k] > 10
#        u[k,gg] = 10
#        
    return e.T, u.T

    
def limitaVelocidade(v, v_lim):
    A = v > v_lim
    v[A] = v_lim
    A = v < -v_lim
    v[A] = -v_lim
    
    
def coefInercial(Wmin, Wmax, i, imax):
    return (Wmax - (Wmax - Wmin)*i/imax)

def atualizaFitness(P,ts,tf,fitpBest,pbest,x,LAMBDA):
    #y = sphereFunction_v2(posAtual)
    #LAMBDA = 1 # coeficiente de supressao de acoes de controle
    #mE, mDU = matrizDeErro(posAtual)
    mE, mDU = picontrol(P,ts,tf,x,len(fitpBest))
    mDU = mDU[...,1:-1] - mDU[...,0:-2]
    f = ise(mE) + LAMBDA*ise(mDU) # MULTIOBJETIVO (LQR)
    #f = ise(mE) # MONO OBJETIVO
    #y2 = ise(mDU)
    A = f < fitpBest # vetor de decisao
    fitpBest[A] = f[A]
    pbest[A] = x[A]
    
def x_bbpso(pBest, fit_pBest, DIM, T_ENXAME):
    
    #np.tile(pBest[np.argmin(fit_pBest)],[T_ENXAME,1]) # Cria uma matriz com GBEST em cada linha
    mean = 0.5*(pBest + np.tile(pBest[np.argmin(fit_pBest)],[T_ENXAME,1]))
    std = np.abs(pBest - np.tile(pBest[np.argmin(fit_pBest)],[T_ENXAME,1]))
    
    return std*np.random.randn(T_ENXAME,DIM) + mean

def atualizaVel(x,v,pbest, gbest, w, c1,c2, k):
    
    T_ENXAME = len(x)
    DIM = len(x[0])
    r1 = np.random.rand(T_ENXAME,DIM) # entre 0 e 1
    r2 = np.random.rand(T_ENXAME,DIM) 
    #w = 0.9
    return k*(w*v + c1*r1*(pbest - x) + c2*r2*(np.tile(gbest,[T_ENXAME,1])-x))


def pso(P,ts,tf,T_ENXAME = 50,    
        Imax = 20, #Numero maximo de iteraçoes para calculo do coef inercial (w)
        Wmax = 0.8,#coef inercial max
        Wmin = 0.4 #coef inercial min   
   ):
    
    
    DIM = 2 #dimensoes do problema    
    
#------ Inicia váriaveis -------------    
    x = 5*np.random.randn(T_ENXAME,DIM)

    v = np.random.randn(T_ENXAME,DIM)    
    pBest = x #pbest recebe o valor de x por ser a única posição conhecida da partic
    fitPbest = np.inf*np.ones(len(x)) 
    fitIter = []
    c1 = 2
    c2 = 2
    #k = 2/np.abs(2-(c1+c2)-np.sqrt((c1+c2)**2-4*(c1+c2)))
    atualizaFitness(P,ts,tf,fitPbest,pBest,x, 0) # atualiza fitness atual e pBest 
    gb = np.argmin(fitPbest)
    fitIter.append(fitPbest[gb])
    for j in range(40):
        
        v = atualizaVel(x,v,pBest,pBest[gb],coefInercial(Wmin,Wmax,j,50),c2,c1,1)  #wIter é um vetor com os valores de w a cada iteraçao
        limitaVelocidade(v, 5000)
        x = x + v
       # limitaRegiao(x, 500)
        atualizaFitness(P,ts,tf,fitPbest,pBest,x,0) # atualiza fitness atual e pBest 
        gb = np.argmin(fitPbest)
        fitIter.append(fitPbest[gb])
          
    return pBest[gb], fitIter
    #return fitIter


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

def resultados(coefs, converg,P,Ts,tf):
    
    t = np.arange(0, tf + Ts, Ts)
    r = 1*np.ones(len(t))


    y1,u1,e1 = picontrol2(P,Ts,tf,coefs[0])
    step_info(t,y1)
    print("ISE: ", np.sum(e1**2))

    fig = plt.figure(figsize=(10,6))
    for i, label in enumerate(['Skogestad IMC', 'IMC', 'Ziegler-Nichols', 'Cohen Coon']):
        ax = fig.add_subplot(2, 2, i+1, ylabel= 'y(t)', xlabel = 't')
        print(label)
        ax.plot(t,y1,label = 'PSO')
        ax.plot(t,r, 'k--')  
        y,u,e = picontrol2(P,Ts,tf,coefs[i+1])
        ax.plot(t,y,label=label)
        step_info(t,y)
        print("ISE: ", np.sum(e**2))

        ax.legend(loc = 'lower right')

    #fig.savefig('respostas.pdf', format='pdf')
   
    fig, ax = plt.subplots(1)
    ax.set_ylabel('ISE')
    ax.set_xlabel('Iterações')
    for i in range(len(converg)):
        ax.plot(converg[i])    
    #fig.savefig('convergencia.pdf', format='pdf')
    
    fig,ax = plt.subplots(1)
    ax.set_ylabel('y(t)')
    ax.set_xlabel('t')
    t,y = scipy.signal.step(P, T=t[0:len(t)//3])
    plt.ylim([0, 1.1])
    ax.plot(t,y, '-')
    #fig.savefig('resposta_degrau.pdf',format = 'pdf')
    step_info(t,y)

    plt.show()

    
def main():
    np.random.seed(0)
    
    DIM = 2 #dimensoes do problema
    T_ENXAME = 50 #tamanho do enxame
    Imax = 50 #Numero maximo de iteraçoes para calculo do coef inercial (w)
    Wmax = 0.7#coef inercial max
    Wmin = 0.2#coef inercial min
    Xlim = 10 # limite do espaço de busca
    
    P = scipy.signal.TransferFunction([1], [1, 4, 6, 4, 1])
   #P = scipy.signal.TransferFunction(np.polymul([1], [-1/2,1]), np.polymul([10,1], [1/2, 1]))

    Ts = 0.1  # time step
    tf = 35 # tempo de simulação
       
    coefs = []
    coefs_aux = []
    fitIter = []
    fitIter_aux = []

    tau = 2.5128
    atraso = 1.8463

    coefs_aux, fitIter_aux = pso(P,Ts,tf)
    coefs.append(coefs_aux)
    fitIter.append(fitIter_aux)
    coefs.append(simc(atraso,tau))
    coefs.append(imc(atraso,tau))
    coefs.append(zn(atraso,tau))
    coefs.append(cc(atraso, tau))
    
    resultados(coefs,fitIter,P,Ts,tf)
 #   for i, lista in enumerate(['pso','simc','imc','zn','cc']):
    print('[Kp, Ki]')
    for i,j in zip(coefs,['pso','simc','imc','zn','cc']):
        print(i,j)


if __name__ == "__main__":
    main()
