#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 21:43:35 2019

@author: Mateus Maruzka Roncaglio
"""
import pickle
import glob
#import pprint
import pypid_control_lib as ctrl

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

def resultados(coefs, converg,P,Ts,tf):
    
    t = np.arange(0, tf + Ts, Ts)
    r = 1*np.ones(len(t))


    fig = plt.figure(figsize = (9,6))
    for i, label in enumerate(['Skogestad IMC', 'IMC', 'Ziegler-Nichols']):
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
#   
#    fig, ax = plt.subplots(1)
#    ax.set_ylabel('ISE')
#    ax.set_xlabel('Iterações')
#    for i in range(len(converg)):
#        ax.plot(converg[i])    
#    #fig.savefig('convergencia.pdf', format='pdf')
#    
#    fig,ax = plt.subplots(1)
#    ax.set_ylabel('y(t)')
#    ax.set_xlabel('t')
#    t,y = scipy.signal.step(P, T=t[0:len(t)//3])
#    plt.ylim([0, 1.1])
#    ax.plot(t,y, '-')
#    #fig.savefig('resposta_degrau.pdf',format = 'pdf')
#    step_info(t,y)

    plt.show()

def plota_simc_imc_zn(fig, ax, y_imc, u_imc, y_simc, u_simc, y_zn, u_zn, t):
    ax[0].step(t,y_imc[0], label = "Internal Model Control")
    ax[0].step(t,y_simc[0], label = "Skogestad IMC")
    ax[0].step(t,y_zn[0], label = "Ziegler-Nichols")
    ax[0].legend()
    
    ax[1].step(t, u_imc[0])
    ax[1].step(t, u_simc[0])
    ax[1].step(t, u_zn[0])
    
def plota_os_pickle(caminho, df, fig, ax):
    i = 0
    for file_name in glob.glob(caminho):
        print(file_name)
        with open(file_name, "rb") as f:
            i = i + 1
            while True:
                try:
                    data = pickle.load(f)
                    params = data.get('Params')
                   
                    y,e,u = ctrl.picontrol(data.get('Process'), params.get('Ts'), params.get('Tf'), np.array([data.get('Gbest')]), 1)
                    t = np.arange(0, params.get('Tf') + params.get('Ts'), params.get('Ts'))

                    os, tr, ts = ctrl.step_info(t,y[0])
                    df.append([data.get('Metodo'),os,tr,ts, ctrl.ise(y), ctrl.tvc(u),params.get('Lambda')])
                   
                    # print("Params: ISE e lambda = %f" %(params.get('Lambda')))
                    # print("Os :%f \nTr: %f\nTs: %f" %(os,tr,ts))

                    r = np.ones(len(t))
                    ax[0].plot(t,r, 'k--', linewidth = 1.2)
                    ax[0].step(t,y[0], label = data.get('Metodo'))
                    ax[0].legend()


                    ax[1].step(t,u[0])

                    ax[0].set_ylabel('y(t)')
                    ax[1].set_ylabel('u(t)')
                    plt.xlabel('t (s)')
                    plt.xlim([0, params.get('Tf')])
                    #plt.show()
                    

#                    df = pd.DataFrame(df)
#                    with open("dados/tabelas/tabelas%i.tex"%i, "w") as g:
#                        g.write(df.to_latex())

                    # plt.plot(np.arange(0, params.get('Ts') + params.get('Tf'), params.get('Ts')),y.T)
                except EOFError:
                    break


def main():
    
    # fopdt = [3.96718829, 1.98359314]
    # fopdt = [1.98359314, 3.96718829]
    fopdt = [17,          2.33357111]
    #fopdt = [1.77033282 ,2.52943189]
    L = 17*0.1
    
    T = 2.33357111
    
    with open("dados/IAE.pickle", "rb") as f:
        data = pickle.load(f)
        
    params = data.get('Params')
    t = np.arange(0, params.get('Tf') + params.get('Ts'), params.get('Ts'))

    df = [["Método","Os", "Tr", "Ts", "ISE", "TVC", "Coeficiente de supressão da ação de controle"]]

    y_imc, e_imc, u_imc = ctrl.picontrol(data.get('Process'), params.get('Ts'), params.get('Tf'), ctrl.imc(L, T), 1)
    y_zn, e_zn, u_zn = ctrl.picontrol(data.get('Process'), params.get('Ts'), params.get('Tf'), ctrl.zn(L, T), 1)
    y_simc, e_simc, u_simc = ctrl.picontrol(data.get('Process'), params.get('Ts'), params.get('Tf'), ctrl.simc(L, T), 1)
       
    os, tr, ts = ctrl.step_info(t,y_imc[0])
    df.append(["IMC",os,tr,ts, ctrl.ise(y_imc), ctrl.tvc(u_imc)])
    os, tr, ts = ctrl.step_info(t,y_zn[0])
    df.append(["Ziegler-Nichols",os,tr,ts, ctrl.ise(y_zn),ctrl.tvc(u_zn)])
    os, tr, ts = ctrl.step_info(t,y_simc[0])
    df.append(["SIMC",os,tr,ts, ctrl.ise(y_simc), ctrl.tvc(u_simc)])
            
    
    fig1, ax = plt.subplots(ncols=1, nrows=2, constrained_layout=True, sharex = True, figsize = (9,6))
    
    plota_os_pickle("dados/IAE*.pickle", df, fig1, ax)
    plota_simc_imc_zn(fig1, ax, y_imc, u_imc, y_simc, u_simc, y_zn, u_zn, t)
    plt.savefig("./graficos/IAE.pdf")

    fig1, ax = plt.subplots(ncols=1, nrows=2, constrained_layout=True, sharex = True, figsize = (9,6))

    plota_os_pickle("dados/ITAE*.pickle", df, fig1, ax)
    plota_simc_imc_zn(fig1, ax, y_imc, u_imc, y_simc, u_simc, y_zn, u_zn, t)
    plt.savefig("./graficos/ITAE.pdf")

    fig1, ax = plt.subplots(ncols=1, nrows=2, constrained_layout=True, sharex = True, figsize = (9,6))

    plota_os_pickle("dados/ISE*.pickle", df, fig1, ax)
    plota_simc_imc_zn(fig1, ax, y_imc, u_imc, y_simc, u_simc, y_zn, u_zn, t)
    plt.savefig("./graficos/ISE.pdf")

    fig1, ax = plt.subplots(ncols=1, nrows=2, constrained_layout=True, sharex = True, figsize = (9,6))
    
    #plota_os_pickle("dados/*.pickle", df, fig1, ax)
    #ax[0].legend()

    #plt.savefig("./graficos/TODAS.pdf")

    i = 0
    df = pd.DataFrame(df)
    with open("dados/tabelas/tabelas%i.tex"%i, "w") as g:
        g.write(df.to_latex())
#    fig2 = plt.figure(constrained_layout=True)
#    spec2 = gridspec.GridSpec(ncols=2, nrows=2, figure=fig2)
#    f2_ax1 = fig2.add_subplot(spec2[0, 0])
#    f2_ax2 = fig2.add_subplot(spec2[0, 1])
#    f2_ax3 = fig2.add_subplot(spec2[1, 0])
#    f2_ax4 = fig2.add_subplot(spec2[1, 1])



    
if __name__ == "__main__":
    main()