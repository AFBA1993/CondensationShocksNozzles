# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 21:53:53 2019

@author: USER
"""
from scipy.optimize import fsolve
from Coollib import Cool_PT, Cool_Ps, Cool_hs, Cool_Psat
from cflow_cool import fnd_isenexp, fnd_throat,pss_Shock, diff_Shockloc, fnd_stagthr
from Geos import Moses_nozzle,Moses_expansion
from Twophase import fnd_2phase_exp, Propmix, sol2Phase
from TwophGeos import Moses_2phase
from Cnd_shock import cshock
from solution import sol_profiles, Deliver_res
from Rankine_line import Rayleighl, Hugoniotl, fnd_CJ, fndPs_RHeq, RayleighlS, fnd_RHS, fnd_CJS
import numpy as np
import matplotlib.pyplot as plt
from momentum import Err_mom
import pandas as pd

"""
Moses 410 To=273.15+104, Po=70727.516 xW=10.020000000000095
Moses 417 To=273.15+106, Po=70020.908 xW=10.2666666666668
Moses 428 To=273.115+100,Po=54702.168 xW=10.340000000000146
Moses 434 To=273.115+100,Po=41356.599 xW=11.225833333333215
"""


f='Water'
EoS='HEOS'

J=1E+20

A2=1
A1=4*A2

To=273.15+100
Po=41356.599
xg=11.225833333333215




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



exp=Moses_expansion(-0.001*P2,A2,P2,s2,f,EoS,xg)
P3=exp[14]
T3=exp[13]
u3=exp[12]

xcond=exp[11]


A3=Moses_nozzle(xcond)*A2

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



cnd=cshock(P3,T3,u3,q,f,EoS)
Psd=cnd[9]
Vsd=1/cnd[10]

Pwd=cnd[5]
Vwd=1/cnd[6]

P4=cnd[5]
X4=cnd[7]
u4=cnd[8]

prop4=Propmix(P4,X4,f,EoS)
s4=prop4[4]



exp_2ph=Moses_2phase(P4,X4,u4,A2,xcond,f,EoS,1000)
sol=sol_profiles(exp,exp_2ph,To,Po)



"""
Folder2send="/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp434/"
Name_Data="434schk"
take=Deliver_res(sol,Folder2send,Name_Data)
"""







