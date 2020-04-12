#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 08:56:31 2019

@author: Mateus Maruzka Roncaglio
"""

import pickle 
import glob
import numpy as np
import matplotlib.pyplot as plt


def main():
    
    
    for file_name in glob.glob("dados/*.pickle"):
        with open(file_name, "rb") as f:
            while True:
                try:
                    data = pickle.load(f)
                    
                    plt.plot(data.get('fitIter'), label = data.get('Metodo'))
                except EOFError:
                    break


 #   plt.ylim([0,1500])
    plt.semilogx()
    plt.semilogy()

    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()