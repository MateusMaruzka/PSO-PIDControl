# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 18:46:58 2020

@author: maruzka
"""


import glob
import pickle

import numpy as np
import matplotlib.pyplot as plt 
import scipy.signal as signal

import pypid_control_lib as pypid



# plt.rcParams["font.family"] = "Times New Roman"

def main():
    
    
    
    P = signal.TransferFunction([1],[1, 2, 4, 2,1 ])
    
    for file_name in glob.glob("dados/*.pickle"):
        
        print(file_name)
        with open(file_name, "rb") as f:
            
            fig, ax = plt.subplots(ncols = 1, nrows = 2, sharex = True)
            
            data = pickle.load(f)
            
            params = data.get("Params")
            
            y,e,u,r,t = pypid.picontrol(P, params.get('Ts'), params.get('Tf'), np.array([data.get('Gbest')]), 1)
             # t = np.arange(0, params.get('Tf') + params.get('Ts'), params.get('Ts'))

            ax[0].step(t,y[0], '--', label = "ISE "+ r'($\lambda$' + ')', linewidth = 1.2)
            ax[0].step(t,r[0], 'k--', label = "Referência", linewidth = 1)
            ax[0].legend()
            
            ax[1].step(t,u[0], '--', linewidth = 1.2)
            
            
            plt.xlim([0, params.get("Tf")])



                
if __name__ == "__main__":
    
    main()