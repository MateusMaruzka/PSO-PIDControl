#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 21:43:35 2019

@author: Mateus Maruzka Roncaglio
"""
import pickle
import glob
import pprint
import pypid_control_lib as ctrl
import matplotlib.pyplot as plt
import numpy as np

from maruzka_plot_lib import maruzka_plot


def main():
    
    for file_name in glob.glob("dados/*.pickle"):
        with open(file_name, "rb") as f:
            while True:
                try:
                    data = pickle.load(f)
                    params = data.get('Params')
                    y,e,u = ctrl.picontrol(data.get('Process'), params.get('Ts'), params.get('Tf'), np.array([data.get('Gbest')]), 1)
                    plt.plot(np.arange(0, params.get('Ts') + params.get('Tf'), params.get('Ts')),y.T)
                    plt.show()
                except EOFError:
                    break
    
    
if __name__ == "__main__":
    main()