# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 13:56:02 2019

@author: sisea
"""
from LibiPRSV import iPRSV_mas_PT, iPRSV_mas_Ps, iPRSV_mas_hs
from cflowiPRSV import fnd_throat, fnd_stagthr, fnd_isenexp 
from Cndshock import cshock
from GeOsiPRSV import Moore_nozzle,Moore_expansion,Moses_expansion, Moses_nozzle
from TphGeos import Moore_2phase, Moses_2phase
from Two_phase import Propmix
from solution import sol_profiles
from Rayleigh_Hugoniot import Hugoniotl, Rayleighl
#from momentum import Err_mom
import pandas as pd
from solution import sol_profiles, Deliver_res
"""
Moses 410 To=273.15+104, Po=70727.516 xW=10.020000000000095
Moses 417 To=273.15+106, Po=70020.908 xW=10.2666666666668
Moses 428 To=273.115+100,Po=54702.168 xW=10.340000000000146
Moses 434 To=273.115+100,Po=41356.599 xW=11.225833333333215
"""


f='Water'
q=0

A2=1
#A1=Moore_nozzle(0)*A2
A1=4*A2

To=273.15+106
Po=70020.908

xg=10.37499999999984

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


exp=Moses_expansion(-0.01*P2,A2,P2,s2,f,xg)
#exp=Moore_expansion(-0.01*P2,A2,P2,s2,f,xg)
P3=exp[14]
T3=exp[13]
u3=exp[12]
xcond=exp[11]


A3=Moses_nozzle(xcond)*A2
#A3=Moore_nozzle(xcond)*A2
prop3=iPRSV_mas_PT(f,P3,T3)
T3=prop3[2]
V3=prop3[1]
h3=prop3[3]
s3=prop3[4]

mass3=A3/V3*u3
e3=h3+0.5*u3**2.0


cnd=cshock(P3,T3,u3,q,f)
Psd=cnd[9]
rhosd=cnd[10]



P4=float(cnd[5])
X4=float(cnd[7])
u4=float(cnd[8])




prop4=Propmix(P4,X4,f)
V4=1/prop4[1]
h4=prop4[3]
s4=prop4[4]

e4=h4+0.5*u4**2.0
mass4=A3*u4/V4

exp2phs=Moses_2phase(P4,X4,u4,1,xcond,f,100)


"""
sol=sol_profiles(exp,exp2phs,To,Po)
Folder2send="/home/andres777/Área de Trabalho/"
Name_Data="iprsvmoses417"
take=Deliver_res(sol,Folder2send,Name_Data)
"""

   





