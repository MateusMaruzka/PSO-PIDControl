#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  4 17:08:41 2019

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


def main():
    
    L = 1.77033282
    T = 2.52943189
    ganhos = [["Método", "Kp", "ki"]]
    
    for file_name in glob.glob("dados/pid_pso_?????*.pickle"):
        with open(file_name, "rb") as f:
            
            while True:
                try:
                    data = pickle.load(f)
                    gbest = data.get('Gbest')
                    ganhos.append([file_name, gbest[0], gbest[1]])

                except EOFError:
                    break
        
    imc = ctrl.imc(L, T);
    ganhos.append(["IMC", imc[0][0], imc[0][1]])
    
    simc = ctrl.simc(L,T)
    ganhos.append(["SIMC", simc[0][0], simc[0][1]])
    
    zn = ctrl.zn(L, T);
    ganhos.append(["Ziegler-Nichols", zn[0][0], zn[0][1]])


    print(ganhos)
    
    ganhos = pd.DataFrame(ganhos)
    with open("dados/tabelas/tabelas_ganhos.tex", "w") as g:
        g.write(ganhos.to_latex())
        
        
if __name__ == "__main__":
    main()
