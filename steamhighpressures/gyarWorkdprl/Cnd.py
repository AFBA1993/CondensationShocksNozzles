# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 20:26:27 2019

@author: USER
"""
import numpy as np
import matplotlib.pyplot as plt
from Coollib import Cool_PT, Cool_Psat, Cool_Tsat
from scipy.optimize import fsolve

f="Water"
EoS="PR"







u1=571.1344128583836
P1=27133.19076511832
T1=292.83477469033056
prop1=Cool_PT(P1,T1,f,EoS)
h1=prop1[3]
rho1=prop1[1]
    
def cndshock(P2,P1,T1,u1,f,EoS):
    prop1=Cool_PT(P1,T1,f,EoS)
    h1=prop1[3]
    rho1=prop1[1]
    propsat=Cool_Psat(P2,f,EoS)

    
#    psat=Cool_Tsat(T2,f,EoS)
#    psat=psat[0] 
    
    rho_g=1/propsat[2]
    h_g=propsat[3]
    

    h_l=propsat[7]
    
    x=(P1-P2+rho1*u1**2.0)*rho_g/(rho1*u1)**2
    y=1-x
      
    h2=x*h_g+y*h_l
    u2=(2*(h1-h2+0.5*u1**2.0))**0.5
    rho2=rho_g/x
    
          
    mass1=rho1*u1
    mass2=rho2*u2
          
    diff=mass1-mass2
    return diff
    
 
    
    
    
    
def correr(f,EoS):    
    u1=571.1344128583836
    P1=27133.19076511832
    T1=292.83477469033056   

    
    
        
    prop1=Cool_PT(P1,T1,f,EoS)
    h1=prop1[3]
    rho1=prop1[1]
    
    
    #Equations for cndshock
    P_trp=612
    Pc=22064000.0-219*1e+05
    Pdiff=Pc-P_trp
    ne=500
    
    err=np.zeros(ne)
    err_e=np.zeros(ne)
    err_m=np.zeros(ne)
    T=np.zeros(ne)
    P2i=np.zeros(ne)
    P=np.zeros(ne)
    xi=np.zeros(ne)
    yi=np.zeros(ne)
    T_rat=np.zeros(ne)
    P_rat=np.zeros(ne)
    
    stp=Pdiff/ne
    P2=P_trp
    
    for i in range (ne):
        
        propsat=Cool_Psat(P2,f,EoS)
        T2=propsat[1]
        
    #    psat=Cool_Tsat(T2,f,EoS)
    #    psat=psat[0] 
        
        rho_g=1/propsat[2]
        h_g=propsat[3]
        
        rho_l=1/propsat[6]
        h_l=propsat[7]
        
        x=(P1-P2+rho1*u1**2.0)*rho_g/(rho1*u1)**2
        y=1-x
          
        h2=x*h_g+y*h_l
        u2=(2*(h1-h2+0.5*u1**2.0))**0.5
        rho2=rho_g/x
        
              
        mass1=rho1*u1
        mass2=rho2*u2
              
        e1=h1+u1**2.0/2
        e2=h2+u2**2.0/2
        mom1=P1+rho1*u1**2.0
        mom2=P2+rho2*u2**2.0
          
    diff=mass1-mass2
    
    
    xi[i]=x
    yi[i]=y
    T[i]=T2
    P[i]=P2/100000
    T_rat[i]=T2/T1
    P_rat[i]=P2/P1
    err[i]=diff
    err_e[i]=e1-e2
    err_m[i]=mom1-mom2
    P2i[i]=P2
    
    P2=P2+stp
    
    return err,err_e,err_m,T_rat, P_rat, xi


#results=correr(f,EoS)

 
    
Pguess=(37537.7+37864.5)/2
Psol=fsolve(cndshock,Pguess,(P1,T1,u1,f,EoS))
tol=cndshock(Psol,P1,T1,u1,f,EoS) 
propsat=Cool_Psat(Psol,f,EoS)
T2=propsat[1]

#    psat=Cool_Tsat(T2,f,EoS)
#    psat=psat[0] 

rho_g=1/propsat[2]
h_g=propsat[3]
c_g=propsat[10]

rho_l=1/propsat[6]
h_l=propsat[7]
c_l=propsat[11]

x=(P1-Psol+rho1*u1**2.0)*rho_g/(rho1*u1)**2
y=1-x
  
h2=x*h_g+y*h_l
u2=(2*(h1-h2+0.5*u1**2.0))**0.5
rho2=rho_g/x

mass1=rho1*u1
mass2=rho2*u2
        
mom1=P1+rho1*u1**2.0
mom2=Psol+rho2*u2**2.0  
e1=h1+u1**2.0/2
e2=h2+u2**2.0/2

v2=1/rho_g*x+y*1/rho_l
c2=x*c_g+y*c_l
M2=u2/c2
rhoxx=1/v2

