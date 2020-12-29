# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 12:14:32 2020

@author: maruzka
"""

import pickle 

f = "dados/ISE/ISE_lteste=0.9900.pickle"


with open(f, "rb") as file:
    
     data = pickle.load(file)
     params = data.get("Params")
     
     print(params)
