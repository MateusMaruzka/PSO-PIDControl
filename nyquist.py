# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 17:04:52 2020

@author: maruzka
"""


import control
import scipy.signal
import numpy as np

from scipy.signal.ltisys import TransferFunction as TransFun
from numpy import polymul,polyadd
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches


# Regular transfer functions from control and scipy.signal don't allow for algebraic operations.
# This class is used to calculate composited transfer functions.
class ltimul(TransFun):
    def tolti(self):
        return scipy.signal.lti(self.num, self.den)

    def __neg__(self):
        return ltimul(-self.num,self.den)

    def __floordiv__(self,other):
        # can't make sense of integer division right now
        return NotImplemented

    def __mul__(self,other):
        if type(other) in [int, float]:
            return ltimul(self.num*other,self.den)
        elif type(other) in [TransFun, ltimul]:
            numer = polymul(self.num,other.num)
            denom = polymul(self.den,other.den)
            return ltimul(numer,denom)

    def __truediv__(self,other):
        if type(other) in [int, float]:
            return ltimul(self.num,self.den*other)
        if type(other) in [TransFun, ltimul]:
            numer = polymul(self.num,other.den)
            denom = polymul(self.den,other.num)
            return ltimul(numer,denom)

    def __rtruediv__(self,other):
        if type(other) in [int, float]:
            return ltimul(other*self.den,self.num)
        if type(other) in [TransFun, ltimul]:
            numer = polymul(self.den,other.num)
            denom = polymul(self.num,other.den)
            return ltimul(numer,denom)

    def __add__(self,other):
        if type(other) in [int, float]:
            return ltimul(polyadd(self.num,self.den*other),self.den)
        if type(other) in [TransFun, type(self)]:
            numer = polyadd(polymul(self.num,other.den),polymul(self.den,other.num))
            denom = polymul(self.den,other.den)
            return ltimul(numer,denom)

    def __sub__(self,other):
        if type(other) in [int, float]:
            return ltimul(polyadd(self.num,-self.den*other),self.den)
        if type(other) in [TransFun, type(self)]:
            numer = polyadd(polymul(self.num,other.den),-polymul(self.den,other.num))
            denom = polymul(self.den,other.den)
            return ltimul(numer,denom)

    def __rsub__(self,other):
        if type(other) in [int, float]:
            return ltimul(polyadd(-self.num,self.den*other),self.den)
        if type(other) in [TransFun, type(self)]:
            numer = polyadd(polymul(other.num,self.den),-polymul(other.den,self.num))
            denom = polymul(self.den,other.den)
            return ltimul(numer,denom)

    # sheer laziness: symmetric behaviour for commutative operators
    __rmul__ = __mul__
    __radd__ = __add__


plt.rcParams["font.family"] = "Times New Roman"

# Process
P = ltimul([1], [1, 4, 6, 4, 1])

# PID Control (Parallel form)
# U(s)/E(s) = Kp + Ki/s + Kd*s
# U(s)/E(s) = Kd*s^2 + Kp*s + Ki)/s

Kp = 1.253
Ki = 0.30815
# Kp=0.941
# Ki=0.273
Kp = 1.167    
Ki = 0.194

Kp = 0.828
Ki = 0.29

Kp = 0.751 
Ki = 0.286

Kd = 0


# Circulo unitario
theta = np.linspace(0, 2*np.pi, 100) 
r = np.sqrt(1.0)
x1 = r*np.cos(theta)
x2 = r*np.sin(theta)



C_gains = {
    "ISE" : [1.253, 0.308],
    "ISE-TVC" : [0.941, 0.273],
    "IAE": [0.924, 0.307],
    "IAE-TVC": [0.828, 0.29],
    "ITAE":[0.751, 0.286],
    "ITAE-TVC":[0.751, 0.286],
    "IMC":[0.648, 0.278],
    "ZN":[1.167, 0.194]
    }


iter_Cgains = iter(C_gains)

dados_tabela = [["Indice", "Margem de Ganho", "Margem de fase", "Sensitivity"]]

# fig = plt.figure(figsize=(20,35))

# widths = [15,15]
# heights = [15,15,15,15]

# gs_kw = dict(width_ratios=widths, height_ratios=heights)
# fig, axes = plt.subplots(ncols=2, nrows=4,
#                           constrained_layout=True,
#                           figsize = (10,30),
#                           sharex = 'all',
#                           sharey = 'all',
#                           gridspec_kw=gs_kw                         )
# fig.tight_layout()

fig, ax = plt.subplots(nrows = 4, ncols = 2, figsize=(3*3.2,12))
ax = ax.flat

# plt.subplots_adjust(left=0.25,
# #                     # bottom=0.2, 
#                     right=0.75, 
# #                     # top=0.8, 
# #                     # wspace=0, 
# #                     # hspace=0)
#                     )

for idx, x in enumerate(C_gains):
    
    _dados = []
    Kp, Ki = C_gains[x]
    Kd = 0
    C = ltimul([Kd, Kp, Ki], [1, 0])
    
    # Closed loop
    G = P*C/(1 + P*C)
    
    # Transform to control library representation for margins
    G2 = control.tf(G.num, G.den)
    G2 = control.minreal(G2)
    print(G2)
    
    print('Closed Loop Roots')
    print(np.roots(G2.den[0][0]))
    
    if np.all( np.roots(np.real(G2.den[0][0])) < 0 ):
        print('Stable system')
    else:
        print('Unstable system') 
    
    gm, pm, sm, _, _, _ = control.stability_margins(G2)
    _dados = [x, gm, pm, 1/sm]
    dados_tabela.append(_dados)
    
    print("Gain \t Phase \t Sensitivity")
    print(gm,'\t', pm,'\t', sm)
    
    real, imag, omega = control.nyquist(G2, linewidth=1, omega = np.arange(0, 2*np.pi, 0.01), Plot = False)
    ax[idx].plot(real, imag, 'b' ,label = x)
    ax[idx].plot(real, -imag,'b' ,label = x)
    ax[idx].set_aspect(1, share = True)

    ax[idx].set_ylabel("Eixo imaginário")
    ax[idx].set_xlabel("Eixo real")
    ax[idx].grid(True)

    # hide ugly arrows
    for arrow in ax[idx].findobj(matplotlib.patches.FancyArrow):
        arrow.set_visible(False) 
    
        
    ax[idx].plot(x1, x2, color = 'r', label = "Círculo\nUnitário")    

    handles, labels = ax[idx].get_legend_handles_labels()
    ax[idx].legend(handles[1:], labels[1:], 
              bbox_to_anchor=(1.05, 1), 
              loc='upper left', 
              borderaxespad=0.,
              fontsize = 9)

    # Limit window size to nyquist circle
    
 
    ax[idx].set_xlim([-1.5, 1.5])
    ax[idx].set_ylim([-1.5, 1.5])
    
    plt.ylabel("Eixo imaginário")
    plt.xlabel("Eixo real")
    # plt.tight_layout()
    

    plt.show()


# right_inset_ax = fig.add_axes([.1, .1, .1, .1], facecolor='k')
# right_inset_ax.set(title='Probability', xticks=[], yticks=[])
# right_inset_ax.set_ylabel('current (nA)')


# plt.ylabel('this is vertical\ntest', multialignment='center')

plt.savefig('nyquist.pdf')
 
df = pd.DataFrame(dados_tabela)
with open("dados/sensitivity.tex", 'w') as tf:
     tf.write(df.to_latex())