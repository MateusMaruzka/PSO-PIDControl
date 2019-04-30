#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 22:01:00 2019

@author: Mateus Maruzka Roncaglio
"""

import matplotlib.pyplot as plt

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
    
    lines = []
    t = iter(colors)
    ax = []
    cols = len(x)
    fig1 = plt.figure(num=1, figsize=figsize)
    for i in y: # caso o eixo x seja o msm para todos
        #line, = plt.step(x,y, '-', color=next(t), linewidth=1.2)
        row = (i // cols)
        col = i % cols
        ax.append(fig1.add_subplot([row, col]))
        ax[-1].plot(x, y, 'o', ls='-')

    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    
    # plt.title(title, loc='center')  
    
    # melhor lugar de legenda
    plt.legend(tuple(lines), tuple(legends), loc = 'best')
    # plt.legend((line6, line1, line2, line3, line4, line5, line6), (legend6, legend1, legend2, legend3, legend4, legend5), loc='best')
    
    plt.savefig("./graficos/"+title+".svg")
    plt.show()



#def resultados(coefs, converg,P,Ts,tf):
#    
#    t = np.arange(0, tf + Ts, Ts)
#    r = 1*np.ones(len(t))
#
#
#    y1,u1,e1 = picontrol2(P,Ts,tf,coefs[0])
#    step_info(t,y1)
#    print("ISE: ", np.sum(e1**2))
#
#    fig = plt.figure(figsize=(10,6))
#    for i, label in enumerate(['Skogestad IMC', 'IMC', 'Ziegler-Nichols', 'Cohen Coon']):
#        ax = fig.add_subplot(2, 2, i+1, ylabel= 'y(t)', xlabel = 't')
#        print(label)
#        ax.plot(t,y1,label = 'PSO')
#        ax.plot(t,r, 'k--')  
#        y,u,e = picontrol2(P,Ts,tf,coefs[i+1])
#        ax.plot(t,y,label=label)
#        step_info(t,y)
#        print("ISE: ", np.sum(e**2))
#
#        ax.legend(loc = 'lower right')
#
#    #fig.savefig('respostas.pdf', format='pdf')
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
#
#    plt.show()

