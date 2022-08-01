# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 22:40:46 2019

@author: USER
"""
from scipy.optimize import fsolve
from Coollib import Cool_PT, Cool_Ps
from cflow_cool import fnd_isenexp, fnd_throat,pss_Shock, diff_Shockloc
from Geos import Moore_nozzle, Barschdorff_nozzle, Barschdorff_expansion, Moore_expansion,Arina_Shock, Arina_expansion



"Initial Conditions of the Problem"
EoS="HEOS"
f='Air'
A2=1
A1=2.5*A2
#A1=Barschdorff_nozzle(-120)*A2
#A1=Moore_nozzle(0)*A2


A5=1.5
T1=288
P1=20*1E+05
P5r=0.8*P1


prop1=Cool_PT(P1,T1,f,EoS)
V1=1/prop1[1]
h1=prop1[3]
s1=prop1[4]
a1=prop1[6]


P2=fnd_throat(0.85*P1,-0.001*P1,A1,A2,V1,h1,s1,f,EoS)
s2=s1
prop2=Cool_Ps(P2,s2,f,EoS)
V2=1/prop2[1]
h2=prop2[3]
a2=prop2[6]
e2=h2+a2**2.0/2
mass2=A2*a2/V2
u1=(2*(e2-h1))**0.5
mass1=A1*u1/V1




r=fsolve(diff_Shockloc,1.25,(P5r,A2,A5,P2,s2,mass2,e2,P2,-0.01*P2,0.7*P1,-0.01*P1,0.005*P2,f,EoS))
A3=r*A2
P3=fnd_isenexp(P2,-0.01*P2,A3,mass2,e2,s2,f,EoS)
s3=s2
Prop3=Cool_Ps(P3,s3,f,EoS)
T3=Prop3[2]
V3=1/Prop3[1]
h3=Prop3[3]
a3=Prop3[6]


u3=(2*(e2-h3))**0.5
M3=u3/a3
e3=h3+u3**2.0/2
mass2=A3*u3/V3

k1=u3/V3
k2=P3+u3**2.0/V3
k3=h3+u3**2.0/2
PT4=pss_Shock(0.7*P1,-0.01*P1,k1,k2,k3,f,EoS)
P4=PT4[0]
T4=PT4[1]

Prop4=Cool_PT(P4,T4,f,EoS)
V4=1/Prop4[1]
h4=Prop4[3]
s4=Prop4[4]
a4=Prop4[6]
u4=k1*V4
M4=u4/a4

mass4=A3*u3/V3
e4=h4+u4**2.0/2


P5=fnd_isenexp(P4,0.05*P4,A5,mass4,e4,s4,f,EoS)
s5=s4

Prop5=Cool_Ps(P5,s5,f,EoS)
V5=1/Prop5[1]
h5=Prop5[3]
a5=Prop5[6]
u5=(2*(e4-h5))**0.5
mass5=A5*u5/V5
e5=h5+u5**2.0/2
M5=u5/a5

#exp=Barschdorff_expansion(-0.02*P2,A2,P2,s2,f,EoS)
#exp=Moore_expansion(-0.02*P2,A2,P2,s2,f,EoS)
exp1=Arina_expansion(-0.02*P2,A2,P2,s2,f,EoS)
exp2=Arina_expansion(+0.02*P2,A2,P2,s2,f,EoS)
shock=Arina_Shock(r,A2,P2,s2,P4,s4,f,EoS)

