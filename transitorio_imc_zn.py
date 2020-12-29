# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 13:52:37 2020

@author: maruzka
"""


import glob
import pickle
import pandas as pd
import numpy as np
import scipy.signal as signal
import pypid_control_lib as pypid


P = signal.TransferFunction([1],[1, 4, 6, 4, 1])





def main():
    
    Kp_zn, Ki_zn = pypid.zn(1.8, 2.3335711, "PI")[0]
    Kp_imc, Ki_imc = pypid.imc(atraso=1.8, tau = 2.3335711, tauc=1.8)[0] # skogestad 2003

    y_zn,e_zn,u_zn,r_zn,t_zn = pypid.picontrol(P, 0.1, 100//2, np.array([[Kp_zn, Ki_zn]]), 1, d=0)
    y_imc,e_imc,u_imc,r_imc,t_imc = pypid.picontrol(P, 0.1, 100//2, np.array([[Kp_imc, Ki_imc]]), 1, d=0)
    
    tabela = []

    tabela.append([Kp_imc, Ki_imc, pypid.step_info(t_imc, y_imc[0])])    
    tabela.append([Kp_zn, Ki_zn, pypid.step_info(t_zn, y_zn[0])])
    
    df = pd.DataFrame(tabela)
    print(df.to_latex())
    


    

if __name__ == "__main__":
    
    main()