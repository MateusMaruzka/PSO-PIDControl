# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 18:46:58 2020

@author: maruzka
"""


import glob
import pickle

import numpy as np

import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker

import scipy.signal as signal

import pypid_control_lib as pypid


plt.rcParams["font.family"] = "Times New Roman"

def main():
    
    P = signal.TransferFunction([1],[1, 4, 6, 4, 1])
    files = glob.glob("dados/*.pickle")
    
    
    fig = plt.figure()
    ax = [fig.add_subplot(221+i) for i in range(len(files))]
    
    for idx, file_name in enumerate(files):
        
        with open(file_name, "rb") as f:
                        
            data = pickle.load(f)
            
            params = data.get("Params")
            
            print(data.get('Gbest'))
            gb= data.get('Gbest')
            y,e,u,r,t = pypid.picontrol(P, params.get('Ts'), params.get('Tf')//2, np.array([gb]), 1, atraso = 0)
            print(np.sum(e[0]**2))
            
            """
            # a,b,c,tp = pypid.step_info(t, y[0])
            # ax[0].annotate(r'$K_p=' + str(round(gb[0],2)) + '$' + '\n'
            #                r'$K_i=' + str(round(gb[1],2)) + '$', fontsize = 8,
            #                xy=(t[tp],y[0][tp]), 
            #                            xytext=(t[tp+10], 2), 
            #                            # arrowprops=dict(facecolor='black', shrink=0.01, width = 0.5, headwidth = 5),
            #                            arrowprops = dict(arrowstyle = '->', 
            #                            connectionstyle="angle3"))
            """
            
            ll = params.get('Lambda')
            label = data.get("Metodo").split("_")[0] + "-TVC" + r' ($\lambda = {:.2f}$)'.format(ll) \
                    if ll != 0 else data.get("Metodo").split("_")[0]
            
            print(label)
            ax[idx].step(t,y[0], '--', label = label, linewidth = 1.2)
            ax[idx].step(t,r[0], 'k--', label = "Referência", linewidth = 1)
            ax[idx].set_ylabel("y(t)")
            ax[idx].set_xlabel("t")
            ax[idx].legend()
            
            # ax[1].step(t,u[0], '--', linewidth = 1.2)
            # ax[1].set_ylabel("u(t)")
            # ax[1].set_xlabel("t")
            
            ax[idx].set_xlim([0, params.get("Tf")//2])
            ax[idx].set_ylim(bottom = 0, top = max(y[0]) + 0.1)
            
            
            # StrMethod formatter
            # setup(ax)
            ax[idx].yaxis.set_major_locator(ticker.MultipleLocator(0.25))
            # ax[idx].yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
            ax[idx].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x}"))
         
            
            
            
            plt.tight_layout()            

                
if __name__ == "__main__":
    
    main()