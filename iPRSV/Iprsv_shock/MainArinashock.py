# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 17:48:28 2019

@author: sisea
"""
import numpy as np
from scipy.optimize import fsolve
from LibiPRSV import iPRSV_mas_PT, iPRSV_mas_Ps,iPRSV_mas_hs 
from cflowiPRSV import  fnd_isenexp, fnd_throat,pss_Shock, diff_Shockloc 
from GeOsiPRSV import Arina_Shock, Arina_expansion

f='Water'
A2=8E-06
A1=2.5*A2
A5=1.5*A2
T1=380
P1=1E+05
P5r=0.75*P1

"Nozzle Principal Points Solutions"
"Nozzle inlet"
prop1=iPRSV_mas_PT(f,P1,T1)
V1=prop1[1]
h1=prop1[3]
s1=prop1[4]
a1=prop1[6]
P2=fnd_throat(0.9*P1,-0.01*P1,A1,A2,V1,h1,s1,f)

"Nozzle Throat"
s2=s1
prop2=iPRSV_mas_Ps(f,P2,s2)
V2=prop2[1]
h2=prop2[3]
a2=prop2[6]
u2=a2
e2=h2+a2**2.0/2
mass2=A2*a2/V2
u1=(2*(e2-h1))**0.5
e1=h1+u1**2.0/2
mass1=A1*u1/V1
M1=u1/a1
M2=u2/a2

"Total conditions at the nozzle inlet and throat"
ho1=h1+u1**2.0/2
propo1=iPRSV_mas_hs(f,ho1,s1,P1)

ho2=h2+u2**2.0/2
propo2=iPRSV_mas_hs(f,ho2,s2,P2)


r=1.32523

#r=fsolve(diff_Shockloc,1.25,(P5r,A2,A5,P2,s2,mass2,e2,f,P2,-0.01*P2,0.7*P1,-0.01*P1,0.005*P2))
A3=A2*r
s3=s2
P3=fnd_isenexp(P2,-0.01*P2,A3,mass2,e2,s2,f)
prop3=iPRSV_mas_Ps(f,P3,s3)
V3=prop3[1]
T3=prop3[2]
h3=prop3[3]
a3=prop3[6]
u3=(2*(e2-h3))**0.5
e3=h3+0.5*u3**2.0
M3=u3/a3
mass3=A3*u3/V3

"Total conditions upstream the normal Shock"
ho3=h3+u3**2.0/2
propo3=iPRSV_mas_hs(f,ho3,s3,P3)

"Downstream the Normal shock"
k1=u3/V3
k2=P3+u3**2.0/V3
k3=h3+u3**2.0/2
PT4=pss_Shock(0.7*P1,-0.01*P1,k1,k2,k3,f)
P4=PT4[0]
T4=PT4[1]
prop4=iPRSV_mas_PT(f,P4,T4)
V4=prop4[1]
h4=prop4[3]
s4=prop4[4]
a4=prop4[6]
u4=k1*V4

M4=u4/a4
e4=h4+u4**2.0/2
mass4=u4*A3/V4

"Total conditions downstream the normal Shock"
ho4=h4+u4**2.0/2
propo4=iPRSV_mas_hs(f,ho4,s4,P4)

"Nozzle outlet"
P5=fnd_isenexp(P4,0.01*P2,A5,mass4,e4,s4,f)
s5=s4
prop5=iPRSV_mas_Ps(f,P5,s5)
V5=prop5[1]
T5=prop5[2]
h5=prop5[3]
a5=prop5[6]
u5=(2*(e4-h5))**0.5

M5=u5/a5
e5=h5+u5**2.0/2
mass5=u5*A5/V5

"Total conditions nozzle outlet"
ho5=h5+u5**2.0/2
propo5=iPRSV_mas_hs(f,ho5,s5,P5)


Schk=Arina_Shock(r,A2,P2,s2,P4,s4,f)
#lol2=Arina_expansion(-0.1*P2,A2,P2,s2,f)