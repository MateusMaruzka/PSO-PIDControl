# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 17:30:10 2020

@author: maruzka
"""

import matplotlib.pyplot as plt

from graphs_and_results import plota

funcs = ('ISE', 'ITAE', 'IAE')

caminhos = ["dados/"+x for x in funcs]


for i,f in zip(caminhos, funcs):
    
    fig = plota(i)
    plt.tight_layout()    
    plt.savefig(i+'/'+ f+'.pdf')

