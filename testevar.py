# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 23:16:07 2020

@author: maruzka
"""

import numpy as np
import matplotlib.pyplot as plt

def main():
    
    variancias = np.zeros(100//5+1)
    
    for i in range(0,100, 5):
        
        x = np.random.chisquare(1,i)
        variancias[i//5] = np.var(x)
        
        
    
    plt.plot(variancias)
    
    
if __name__ == "__main__":
    main()