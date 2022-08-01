# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 21:53:53 2019

@author: USER
"""
from scipy.optimize import fsolve
from Coollib import Cool_PT, Cool_Ps, Cool_hs, Cool_Psat, Cool_Tsat
from cflow_cool import fnd_isenexp, fnd_throat_MP, fnd_stagthr,fnd_throat
from Geos import BierB2_expansion, B2_Bier
from drpl_growth import Gyarmathy_growth 
from solution import sol_profiles, Deliver_res
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

f='CarbonDioxide'
EoS='HEOS'


k=1.00
J=1E+20
xg=1.0

A2=1
A1=4*A2



To=300.12
Po=45*100000


q=0


propo=Cool_PT(Po,To,f,EoS)
ho=propo[3]
so=propo[4]


P1=fnd_stagthr(0.995*Po,A1,A2,ho,so,f,EoS,To)
#P1=5818608.44166945
s1=so


prop1=Cool_Ps(P1,s1,f,EoS,To)
V1=1/prop1[1]
T1=prop1[2]
h1=prop1[3]
s1=prop1[4]


P2=fnd_throat_MP(P1,A1,A2,V1,h1,s1,f,EoS,To)

s2=s1

prop2=Cool_Ps(P2,s2,f,EoS,To)

V2=1/prop2[1]
T2=prop2[2]
h2=prop2[3]
s2=prop2[4]
u2=prop2[6]

u1=(2*(h2+0.5*u2**2-h1))**0.5
Psat=Cool_Tsat(T2,f,EoS)
S2=P2/Psat[0]


mass1=A1/V1*u1
mass2=A2/V2*u2

e1=h1+0.5*u1**2.0
e2=h2+0.5*u2**2.0

Propo2=Cool_hs(e2,s2,Po,f,EoS,To+5)


exp=BierB2_expansion(-0.001*P2,A2,P2,s2,f,EoS,J,To,xg,k)


P3=exp[14]
T3=exp[13]
u3=exp[12]

xcond=exp[11]
A3=B2_Bier(xcond)*A2



prop3=Cool_PT(P3,T3,f,EoS)
V3=1/prop3[1]
h3=prop3[3]
s3=prop3[4]
#propcrit=Cool_PT(7377300.0,304.1282,f,EoS)
#scrit=propcrit[4]



mass3=A3/V3*u3
e3=h3+0.5*u3**2.0

expdprl=Gyarmathy_growth(P3,T3,u3,A2,10000,xcond,f,EoS,k)
sol=sol_profiles(exp, expdprl, To, Po)



Folder2send="/media/andres/DISPOSITIVO/Mestrado/Andressss/final2/Bierco2/mesh/10000/"
Name_Data="bierdrpl10000"
Deliver_res(sol,Folder2send,Name_Data)


"""
l=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/final2/RESULTADOSFINAL/Simulaciones/bierco2/bierdrpl4000_tl.csv")
sk=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/final2/Bierco2results/bierco2schk/co2shckbierco2schk_schk.csv")

plt.plot(l["xTd"],l['yTd'],"--")
plt.plot(sk["x"],sk['T'],"--")
"""


