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
import matplotlib
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


def main():
    
    fopdt = [3.96718829, 1.98359314]
    i = 0
    for file_name in glob.glob("dados/*.pickle"):
        with open(file_name, "rb") as f:
            df = [["Os", "Tr", "Ts"]]
            i = i + 1
            while True:
                try:
                    data = pickle.load(f)
                    params = data.get('Params')
                    t = np.arange(0, params.get('Tf') + params.get('Ts'), params.get('Ts'))
                    y,e,u = ctrl.picontrol(data.get('Process'), params.get('Ts'), params.get('Tf'), np.array([data.get('Gbest')]), 1)
                    os, tr, ts = ctrl.step_info(t,y[0])
                    df.append([os,tr,ts])

                    print("Params: ISE e lambda = %f" %(params.get('Lambda')))
                    print("Os :%f \nTr: %f\nTs: %f" %(os,tr,ts))
                    y_imc, e_imc, u_imc = ctrl.picontrol(data.get('Process'), params.get('Ts'), params.get('Tf'), ctrl.imc(fopdt[1], fopdt[0]), 1)
                    y_zn, e_zn, u_zn = ctrl.picontrol(data.get('Process'), params.get('Ts'), params.get('Tf'), ctrl.zn(fopdt[1], fopdt[0]), 1)
                    y_simc, e_simc, u_simc = ctrl.picontrol(data.get('Process'), params.get('Ts'), params.get('Tf'), ctrl.simc(fopdt[1], fopdt[0]), 1)
                    os, tr, ts = ctrl.step_info(t,y_zn[0])
                    df.append([os,tr,ts])

                    fig1, ax = plt.subplots(ncols=1, nrows=2, constrained_layout=True, sharex = True)
                    
                    
                    r = np.ones(len(t))
                    ax[0].plot(t,r, 'k--', linewidth = 1.2)
                    ax[0].plot(t,y[0], label = "PSO")
                    ax[0].plot(t,y_imc[0], label = "IMC")
                    ax[0].plot(t,y_zn[0] , label = "Ziegler-Nichols")
                    ax[0].plot(t,y_simc[0], label = "SIMC")
                    ax[0].legend()


                    ax[1].plot(t,u[0])
                    ax[1].plot(t,u_imc[0])
                    ax[1].plot(t,u_zn[0])
                    ax[1].plot(t,u_simc[0])


                    ax[0].set_ylabel('y(t)')
                    ax[1].set_ylabel('u(t)')
                    plt.xlabel('t (s)')
                    plt.xlim([0, params.get('Tf')])
                    plt.savefig("./graficos/%i.svg" %(i))
                    plt.show()
                    

                    df = pd.DataFrame(df)
                    with open("dados/tabelas/tabelas%i.tex"%i, "w") as g:
                        g.write(df.to_latex())

                    # plt.plot(np.arange(0, params.get('Ts') + params.get('Tf'), params.get('Ts')),y.T)
                except EOFError:
                    break
    

    
#    fig2 = plt.figure(constrained_layout=True)
#    spec2 = gridspec.GridSpec(ncols=2, nrows=2, figure=fig2)
#    f2_ax1 = fig2.add_subplot(spec2[0, 0])
#    f2_ax2 = fig2.add_subplot(spec2[0, 1])
#    f2_ax3 = fig2.add_subplot(spec2[1, 0])
#    f2_ax4 = fig2.add_subplot(spec2[1, 1])



    
if __name__ == "__main__":
    main()