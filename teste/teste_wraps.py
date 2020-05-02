# -*- coding: utf-8 -*-
"""
Created on Fri May  1 19:02:27 2020

@author: maruzka
"""

from functools import wraps



def hehehe():
    
    print("hehehe")
    
    
@wraps(hehehe)
def main():
    
    print("main")   
    
    
    