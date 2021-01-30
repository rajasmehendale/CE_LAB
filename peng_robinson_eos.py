# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 15:20:26 2021

@author: Rajas
"""

import numpy as np
import scipy.optimize as sco
import math
import pandas as pd


class peng_robinson:
        def __init__(self):
            self.tc = 369.3 #kelvin
            self.pc = 49.1*1e5 #bar
            self.w = 0.221
            self.t = 373
            self.r = 8.314472
            self.p_guess = 101325
        
        def volumes(self,p):
            t= self.t
            tc=self.tc
            pc=self.pc
            w=self.w
            r=self.r
            k = 0.37464+(1.54226*w)-(0.26992*pow(w,2))
            root_alpha = 1+(k*(1-math.sqrt(t/tc)))
            self.alpha=pow(root_alpha,2)
            self.a=((0.45724*pow((r*tc),2))/pc)*self.alpha
            self.b=(0.07780*r*tc)/pc    
            a = self.a
            r=self.r
            t=self.t
            b=self.b
            def eqn_solved(v):
                v1=v[0]
                v2=v[1]
                F=np.empty((2))
                F[0] = p - ((r*t)/(v1-b)) + (a/( v1*(v1+b) + b*(v1-b) ))
                F[1] = p - ((r*t)/(v2-b)) + (a/( v2*(v2+b) + b*(v2-b) ))
                return F          
            self.guess = np.array([(r*t/p),1.2*b])
            #v1= sco.fsolve(eqn_solved,guess[0])
            #print(v1)
            #v2= sco.fsolve(eqn_solved,guess[1])
            self.vol= sco.fsolve(eqn_solved,self.guess)
            vg = self.vol[0]
            vl = self.vol[1]
            self.gres_g = (p*vg-r*t) - r*t*np.log(p*(vg-b)/(r*t)) +(a/(b*2*np.sqrt(2)))*np.log((vg+(1-np.sqrt(2))*b)/(vg+(1+np.sqrt(2))*b))
        
            self.gres_l = (p*vl-r*t) - r*t*np.log(p*(vl-b)/(r*t)) +(a/(b*2*np.sqrt(2)))*np.log((vl+(1-np.sqrt(2))*b)/(vl+(1+np.sqrt(2))*b))
            error = self.gres_g  - self.gres_l
            
            self.hres_g = (p*vg-r*t) +((a*(1-np.sqrt(t/tc)))/(b*2*np.sqrt(2)))*np.log((vg+(1-np.sqrt(2))*b)/(vg+(1+np.sqrt(2))*b))
            self.hres_l = (p*vl-r*t) +((a*(1-np.sqrt(t/tc)))/(b*2*np.sqrt(2)))*np.log((vl+(1-np.sqrt(2))*b)/(vl+(1+np.sqrt(2))*b))
        
            return error
##################################################################  
aa = peng_robinson()
aa.tc = 647.14
aa.pc = 220.64*10**5
aa.w = 0.344
#aa.t=373
################################################################## 
t_list = [80.0,85.0,90.0,95.0,99.0]
print("\n ----------------------------\n")
for i in t_list:
    err=1
    aa.t  =i+273.16
    p1=101325.0*0.5
    p2=101325.0*1.0
    err1 = aa.volumes(p1)
    #print("P = %.3f bar, Gres_L = %3.2f, Gres_G = %3.2f" %(p1/1e5,aa.gres_l,aa.gres_g)) 
    err2 = aa.volumes(p2)
    #print("P = %.3f bar, Gres_L = %3.2f, Gres_G = %3.2f" %(p2/1e5,aa.gres_l,aa.gres_g)) 
    p_new= ((p1-p2)/(err1-err2))*(-err1) + p1
    while abs(err)>0.001:
        err =aa.volumes(p_new)
        p1=p2
        p2=p_new
        err1= err2
        err2= err
        p_new= ((p1-p2)/(err1-err2))*(-err1) + p1
    print("Temperature = %.2f K\n" %(aa.t)) 
    data = [['P_Bar', '%.4f'%(p_new/1e5)],['G_res', '%.2f'%(aa.gres_l)],['H_res_L', '%.2f'%(aa.hres_l)],['H_res_G', '%.3f'%(aa.hres_g)]  ]
    df = pd.DataFrame(data, columns = ['Type','Value'])
    display(df)
    print("\n ----------------------------\n")
    
# err=1
# p1=101325.0*0.5
# p2=101325.0*1
# err1 = aa.volumes(p1)
# #print("P = %.3f bar, Gres_L = %3.2f, Gres_G = %3.2f" %(p1/1e5,aa.gres_l,aa.gres_g)) 
# err2 = aa.volumes(p2)
# #print("P = %.3f bar, Gres_L = %3.2f, Gres_G = %3.2f" %(p2/1e5,aa.gres_l,aa.gres_g)) 
# p_new= ((p1-p2)/(err1-err2))*(-err1) + p1
# while abs(err)>0.001:
#     err =aa.volumes(p_new)
#     p1=p2
#     p2=p_new
#     err1= err2
#     err2= err
#     p_new= ((p1-p2)/(err1-err2))*(-err1) + p1
# print("Temperature = %.2f K\n" %(aa.t)) 
# data = [['P_Bar', '%.4f'%(p_new/1e5)],['G_res', '%.2f'%(aa.gres_l)],['H_res_L', '%.2f'%(aa.hres_l)],['H_res_G', '%.3f'%(aa.hres_g)]  ]
# df = pd.DataFrame(data, columns = ['Type','Value'])
# display(df)
# print("\n ----------------------------\n")
