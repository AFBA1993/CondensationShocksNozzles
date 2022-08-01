# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 21:53:53 2019

@author: USER
"""
from scipy.optimize import fsolve
from Coollib import Cool_PT, Cool_Ps, Cool_hs, Cool_Psat, Cool_Tsat
from cflow_cool import fnd_isenexp, fnd_throat, fnd_stagthr
from Geos import B2M_expansion, B2M,B4M_expansion, B4M
from Twophase import fnd_2phase_exp, Propmix, sol2Phase
from TwophGeos import Barschdorff_2phase, B2M_2phase, B4M_2phase
from Cnd_shock import cshock
from solution import sol_profiles, Deliver_res
from Rankine_line import Rayleighl, Hugoniotl, fnd_CJ, fndPs_RHeq, RayleighlS, fnd_RHS, fnd_CJS
import numpy as np
import matplotlib.pyplot as plt
from momentum import Err_mom
import pandas as pd
"""
Total conditions

B4M A1=4.2275*A2
B4M-18C Po=100.7 To=342.2, IF95lim=400, P2=0.6322*Po xg=10
B4M-18B Po=100.7 To=365.53, IF95lim=400

"""


f='Water'
EoS='HEOS'



J=1E+20

A2=1
A1=4.2275*A2

xg=9.320762000000128

To=273.15+365.53
Po=100.7*100000

lim=400
q=0

propo=Cool_PT(Po,To,f,EoS)
ho=propo[3]
so=propo[4]

P1=fnd_stagthr(Po,A1,A2,ho,so,f,EoS,To,lim)
P1=float(P1)
s1=so


prop1=Cool_Ps(P1,s1,f,EoS,To,lim)

V1=1/prop1[1]
T1=prop1[2]
h1=prop1[3]
s1=prop1[4]

P2=0.6322*Po

#P2=6023772.88448448
#s2=s1
s2=s1
prop2=Cool_Ps(P2,s2,f,EoS,To,lim)
T2=prop2[2]
Psat2=Cool_Tsat(T2,f,EoS)
Psat2=Psat2[0]
Subc2=P2/Psat2
print (Subc2)



prop2=Cool_Ps(P2,s2,f,EoS,To,lim)

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


exp=B4M_expansion(-0.001*P2,A2,P2,s2,f,EoS,J,To,lim,xg)



P3=exp[14]
T3=exp[13]
u3=exp[12]

xcond=exp[11]
A3=B4M(xcond)*A2





prop3=Cool_PT(P3,T3,f,EoS)
V3=1/prop3[1]
h3=prop3[3]
s3=prop3[4]

mass3=A3/V3*u3
e3=h3+0.5*u3**2.0



cndRHS=fnd_RHS(P3,T3,u3,1,f,EoS)
Vwd1=cndRHS[0]
Pwd1=cndRHS[1]

Vsd1=cndRHS[3]
Psd1=cndRHS[4]


Ps=fndPs_RHeq(P3,T3,u3,q,f,EoS)


cnd=cshock(P3,T3,u3,q,f,EoS)

#Psd=cnd[9]
#Vsd=1/cnd[10]

Pwd=cnd[5]
Vwd=1/cnd[6]

P4=cnd[5]
X4=cnd[7]
u4=cnd[8]

prop4=Propmix(P4,X4,f,EoS)
s4=prop4[4]


exp_2ph=B4M_2phase(P4,X4,u4,A2,xcond,f,EoS)
sol=sol_profiles(exp,exp_2ph,To,Po)

#plt.plot(sol[0],sol[11],"-")

"""
Folder2send="/home/andres/Área de Trabalho/"
Name_Data="Gyarmth18bshck"
take=Deliver_res(sol,Folder2send,Name_Data)
"""

