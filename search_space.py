# -*- coding: utf-8 -*-
"""
Created on Fri May  1 21:57:48 2020

@author: maruzka
"""


import scipy.signal as signal
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.colors as colors


import pypid_control_lib as pypid




def main():
    
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = 'stix'

    P = signal.TransferFunction([1],[1, 4, 6, 4, 1])
    Ts = 0.1
    Tf = 100
    
    kp = np.arange(-0.1, 4.5, 0.025)
    ki = np.arange(-0.1, 4.5, 0.025)
    
    KP, KI = np.meshgrid(kp,ki)
    # print(KP)
    # print(KI)
    
    Z = []
    
    fig = plt.figure(figsize=(9,9))
    ax = fig.add_subplot(211, projection='3d')

    
    for kps, kis in zip(KP,KI):
        _, e, u, _, _, =  pypid.picontrol(P, Ts, Tf, np.array([[kps,kis]]), num_controladores=len(kp))
        # aux = pypid.ise(e)
        # aux[aux > 900] = 900
        # ax.scatter(kps, kis, aux)
        # print(kps, kis, pypid.ise(e))
        Z.append(pypid.ise(e))
            
        
        
        
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    
    Z = np.array(Z)
    # print(Z)
    # Z[Z > 900] = 900
    
  
    ax.plot_surface(KP, KI, np.log10(Z), 
                    norm = colors.LogNorm(vmin = np.log10(Z.min()), 
                                          vmax = np.log10(Z.max())),
                    cmap='RdBu_r')
    ax.set_xlabel("$K_p$")
    ax.set_ylabel("$K_i$")
    ax.set_zlabel(r'$\log_{10}(ISE)$')
    
    
    
    
    ax2 = fig.add_subplot(212)
    ax2.set_xlabel("$K_p$")
    ax2.set_ylabel("$K_i$")
   
    # ax2.contourf(KP, KI, Z, cmap = 'viridis')
    
    # Z[Z > 900] = 900
    pcm = ax2.pcolor(KP, KI, np.log10(Z),
                    norm=colors.LogNorm(vmin = np.log10(Z.min()), 
                                        vmax= np.log10(Z.max())),
                    cmap='RdBu_r')
    
    
    vmin = Z.min()
    kp_min, ki_min = np.where(Z == vmin)
    ax2.scatter(kp[kp_min], ki[ki_min], np.log10(Z[kp_min, ki_min]))
    
    ax2.set_xlim([0, 4.5])
    ax2.set_ylim([0, 4.5])
    
    cb =  fig.colorbar(pcm, ax=ax2, extend='max')
    
    cb.set_label(r'$\log_{10}(ISE)$')

    plt.tight_layout()
    plt.show()
    
    



if __name__ == "__main__":
    main()