# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 14:30:50 2020

@author: maruzka
"""
import numpy as np
import scipy.signal 
import control.matlab
import control
import matplotlib.pyplot as plt
import matplotlib.patches


from scipy.signal.ltisys import TransferFunction as TransFun
from numpy import polymul,polyadd


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

def main():
    
    Ts = 0.1
    P = scipy.signal.TransferFunction([1], [1, 4, 6, 4, 1])
    Pd = P.to_discrete(Ts);
    w, H = scipy.signal.dfreqresp(Pd)

    Kp = 0.924
    Ki = 0.307
    Cpid = scipy.signal.TransferFunction([(Kp+Ki*Ts), -Kp],[1,-1], dt = Ts)
    # Cpid = scipy.signal.TransferFunction([-Kp, (Kp+Ki*Ts)],[-1,1], dt = Ts)

    # fig, ax = plt.subplots(figsize=(5,5))
    # circ = plt.Circle((0, 0), radius=1, edgecolor='black', facecolor='None', fill = False)
    # ax.plot(H.real, H.imag, "b")
    # ax.plot(H.real, -H.imag, "r")
    # ax.add_artist(circ)
    
    # ax.scatter(-1, 0, marker = 'o')
   
    # ax.set_xlim([-1.1,1.1])
    # ax.set_ylim([-1.1,1.1])
    
    
    # Process
    P = ltimul(Pd.num, Pd.den)
   
    # PID Control (Parallel form)
    # U(s)/E(s) = Kp + Ki/s + Kd*s
    # U(s)/E(s) = Kd*s^2 + Kp*s + Ki)/s

    C = ltimul(Cpid.num, Cpid.den)
    
    # Closed loop
    # G = P*C/(1 + P*C)
    G = 1 + P*C
    # G = P
    
    # Transform to control library representation for margins
    G2 = control.tf(G.num, G.den, Ts)
    G2 = control.minreal(G2)
    print(G2)
    
    print('Closed Loop Roots')
    print(np.roots(G2.den[0][0]))
        

    if np.all( np.roots(np.real(G2.den[0][0])) < 0 ):
        print('Stable system')
    else:
        print('Unstable system') 
    
    gm, pm, sm, _, _, _ = control.stability_margins(G2)
    
    print("Gain \t Phase \t Sensitivity")
    print(gm,'\t', pm,'\t', sm)
    
    # Nyquist plot
    fig, ax = plt.subplots()
    control.nyquist_plot(G2, linewidth=1)
    
    # roots = np.roots(G2.den[0][0])
    # print(roots)
    # [ax.scatter(x.real, x.imag) for x in roots]

    
    # hide ugly arrows
    for arrow in ax.findobj(matplotlib.patches.FancyArrow):
        arrow.set_visible(False)
    
    # Draw a red circle from origin to (-1,1) (Nyquist point at -1)
    nyquistCircle = plt.Circle((0, 0), 1, color='r', fill=False)
    ax.add_artist(nyquistCircle)
    
    # Limit window size to nyquist circle
    ax.set_xlim((-1, 1))
    ax.set_ylim((-1, 1))
    
    plt.show()
        


if __name__ == "__main__":
    main()
