# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 13:56:02 2019

@author: sisea
"""
from LibiPRSV import iPRSV_mas_PT, iPRSV_mas_Ps, iPRSV_mas_hs
from cflowiPRSV import fnd_throat, fnd_stagthr, fnd_isenexp 
from GeOsiPRSV import Moore_nozzle,Moore_expansion, Moses_expansion, Moses_nozzle
from drpl_growth import Chinos_growth, Gyarmathy_growth_Moses ,Gyarmathy_growth_Moore
import pandas as pd
import time
from solution import sol_profiles, Deliver_res


"""
Moore To=358.11, Po=25000
Moses 417 To=273.15+106, Po=70020.908 k=0.62 
 
"""


f='Water'
k=0.62
xg=8.05


A2=1

#A1=Moore_nozzle(0)*A2
A1=4*A2

To=273.15+106
#Po=70727.516
Po=70020.908 

propo=iPRSV_mas_PT(f,Po,To)
ho=propo[3]
so=propo[4]


P1=fnd_stagthr(Po,A1,A2,ho,so,f)
s1=so

prop1=iPRSV_mas_Ps(f,P1,s1)
T1=prop1[2]
V1=prop1[1]
h1=prop1[3]




P2=fnd_throat(P1,-0.005*P1,A1,A2,V1,h1,s1,f)
s2=s1
prop2=iPRSV_mas_Ps(f,P2,s2)
T2=prop2[2]
V2=prop2[1]
h2=prop2[3]
u2=prop2[6]


u1=(2*(h2+0.5*u2**2-h1))**0.5

mass1=A1/V1*u1
mass2=A2/V2*u2

e1=h1+0.5*u1**2.0
e2=h2+0.5*u2**2.0



#exp=Moore_expansion(-0.01*P2,500,A2,P2,s2,f,k,xg)
exp=Moses_expansion(-0.01*P2,500,A2,P2,s2,f,k,xg)

P3=exp[14]
T3=exp[13]
u3=exp[12]
J3=exp[17]
rcrit3=exp[18]
xcond=exp[11]

A3=Moses_nozzle(xcond)*A2

prop3=iPRSV_mas_PT(f,P3,T3)
T3=prop3[2]
V3=prop3[1]
h3=prop3[3]
s3=prop3[4]

mass3=A3/V3*u3
e3=h3+0.5*u3**2.0

expdrpl=Gyarmathy_growth_Moses(P3,T3,u3,A2,4000,xcond,f,k)
sol=sol_profiles(exp,expdrpl,To,Po)


Folder2send="/media/andres/DISPOSITIVO/Mestrado/Andressss/final2/RESULTADOSFINAL/Simulaciones/bierco2"
Name_Data="mosesdprliprsv4000"
take=Deliver_res(sol,Folder2send,Name_Data)   




