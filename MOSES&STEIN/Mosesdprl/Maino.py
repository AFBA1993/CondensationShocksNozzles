# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 21:53:53 2019

@author: USER
"""
from scipy.optimize import fsolve
from Coollib import Cool_PT, Cool_Ps, Cool_hs, Cool_Psat, Cool_DT
from cflow_cool import fnd_isenexp, fnd_throat,pss_Shock, diff_Shockloc, fnd_stagthr
from Geos import Moses_nozzle, Moses_expansion
from drpl_growth import Gyarmathy_growth
from solution import sol_profiles, Deliver_res
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



"""
Moses 410 To=273.15+104, Po=70727.516 k=0.67 
Moses 417 To=273.15+106, Po=70020.908 k=0.665
Moses 428 To=273.115+100,Po=54702.168 k=0.61
Moses 434 To=273.115+100,Po=41356.599 k=0.545
"""


f='Water'
EoS='HEOS'

xg=16
k=0.545
J=1E+20

A2=1


#A1=Moore_nozzle(0)*A2
A1=4*A2

To=273.115+100
Po=41356.599

q=0

propo=Cool_PT(Po,To,f,EoS)
ho=propo[3]
so=propo[4]

P1=fnd_stagthr(Po,A1,A2,ho,so,f,EoS)
s1=so

prop1=Cool_Ps(P1,s1,f,EoS)

V1=1/prop1[1]
T1=prop1[2]
h1=prop1[3]
s1=prop1[4]

P2=fnd_throat(P1,-0.005*P1,A1,A2,V1,h1,s1,f,EoS)
s2=s1


prop2=Cool_Ps(P2,s2,f,EoS)

V2=1/prop2[1]
T2=prop2[2]
h2=prop2[3]
s2=prop2[4]
u2=prop2[6]

u1=(2*(h2+0.5*u2**2-h1))**0.5

mass1=A1/V1*u1
mass2=A2/V2*u2

e1=h1+0.5*u1**2.0
e2=h2+0.5*u2**2.0
exp=Moses_expansion(-0.005*P2,A2,P2,s2,f,EoS,J,xg,k)

P3=exp[14]
T3=exp[13]
u3=exp[12]
J3=exp[17]
rcrit3=exp[18]
xcond=exp[11]

prop3=Cool_PT(P3,T3,f,EoS)
V3=1/prop3[1]

A3=Moses_nozzle(xcond)*A2
mass3=u3*A3/u3

expdrpl=Gyarmathy_growth(P3,T3,u3,A2,1000,xcond,f,EoS,k)
sol=sol_profiles(exp,expdrpl,To,Po)



"""
Folder2send="/media/andres777/DISPOSITIVO/Mestrado/Andressss/final2/RESULTADOSFINAL/Simulaciones/resmoses/exp434/"
Name_Data="Gyar434drpl"
take=Deliver_res(sol,Folder2send,Name_Data)
"""


"""
g1=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/final2/RESULTADOSFINAL/Simulaciones/resmoses/exp410/Gyar410drpl_g1.csv")
g2=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/final2/RESULTADOSFINAL/Simulaciones/resmoses/exp410/Gyar410drpl_g2.csv")
rm=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/final2/RESULTADOSFINAL/Simulaciones/resmoses/exp410/Gyar410drpl_rm.csv")
l=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/final2/RESULTADOSFINAL/Simulaciones/resmoses/exp410/Gyar410drpl_tl.csv")


plt.plot(l["xTd"],l['yTd'],"--")


#plt.plot(g1["xg1"],g1['yg1'],"--")
#plt.plot(g2["xg2"],g2['yg2'],"--")
#plt.plot(rm["xnew"],rm['ynew'])
"""







