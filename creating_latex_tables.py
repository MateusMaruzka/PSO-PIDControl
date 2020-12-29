# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 18:17:42 2020

@author: maruzka
"""

import glob
import pickle
import pandas as pd
import numpy as np
import scipy.signal as signal
import pypid_control_lib as pypid


metodo = "IAE"

files = glob.glob("dados" + "/" +metodo+"/*.pickle")


P = signal.TransferFunction([1],[1, 4, 6, 4, 1])
    
# Kp_zn, Ki_zn = pypid.zn(1.8, 2.3335711, "PI")[0]
# Kp_imc, Ki_imc = pypid.imc(atraso=1.8, tau = 2.3335711, tauc=1.8)[0] # skogestad 2003
# y_zn,e_zn,u_zn,r_zn,t_zn = pypid.picontrol(P, 0.1, 100//2, np.array([[Kp_zn, Ki_zn]]), 1, d=0.0001)
# y_imc,e_imc,u_imc,r_imc,t_imc = pypid.picontrol(P, 0.1, 100//2, np.array([[Kp_imc, Ki_imc]]), 1, d=0.0001)


func = {
     "ISE" : pypid.ise,
     "ITAE" : pypid.itae,
     "IAE" : pypid.iae,
     "TVC" : pypid.tvc
     }

dados = []

for idx, file_name in enumerate(files):
    
    with open(file_name, "rb") as f:
                    
        data = pickle.load(f)
        
        params = data.get("Params")
        gb = data.get('Gbest')
        y, e, u, r, t = pypid.picontrol(P=P, ts=params.get("Ts"),
                                        tf=params.get('Tf'),
                                        vetor_ganhos=np.array([gb]), num_controladores=1)
        
        os, tr, ts = pypid.step_info(t, y[0])
        
        ll = params.get('Lambda')
        label = data.get("Metodo").split("_")[0] + "-TVC" + r' ($\lambda = {:.4f}$)'.format(ll) \
                    if ll != 0 else data.get("Metodo").split("_")[0]
                    
        __dados = []
        __dados.append(label)
        __dados.append(gb)
        method = func[data.get("Metodo").split("_")[0]]
        m = func["TVC"]
        __dados.append(method(e))
        __dados.append(m(u))
        __dados.append(os)
        __dados.append(tr)
        __dados.append(ts)
        __dados.append(ll)
        
        dados.append(__dados)
        
        
       



Kp_zn, Ki_zn = pypid.zn(1.8, 2.3335711, "PI")[0]
Kp_imc, Ki_imc = pypid.imc(atraso=1.8, tau = 2.3335711, tauc=1.8)[0] # skogestad 2003

y_zn,e_zn,u_zn,r_zn,t_zn = pypid.picontrol(P, 0.1, 100//2, np.array([[Kp_zn, Ki_zn]]), 1, d=0)
y_imc,e_imc,u_imc,r_imc,t_imc = pypid.picontrol(P, 0.1, 100//2, np.array([[Kp_imc, Ki_imc]]), 1, d=0)

tabela = []

dados.append(['IMC', Kp_imc, Ki_imc, func[metodo](e_imc), func["TVC"](u_imc), pypid.step_info(t_imc, y_imc[0])])    
dados.append(["ZN",Kp_zn, Ki_zn, func[metodo](e_zn), func["TVC"](u_zn), pypid.step_info(t_zn, y_zn[0])])

df = pd.DataFrame(dados)
with open('dados/' +metodo+ "/" + metodo + ".tex", 'w') as tf:
     tf.write(df.to_latex())
        