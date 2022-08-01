# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 20:26:27 2019

@author: USER
"""
import numpy as np
import matplotlib.pyplot as plt
from Coollib import Cool_PT, Cool_Psat, Cool_Tsat

def errshock(P1,T1,u1,f,EoS):
    
    prop1=Cool_PT(P1,T1,f,EoS)
    h1=prop1[3]
    rho1=prop1[1]
    
    
    #Equations for cndshock
    P_trp=612
    Pc=22064000.0-100000
    Pdiff=Pc-P_trp
    ne=500
    
    err=np.zeros(ne)
    T=np.zeros(ne)
    P=np.zeros(ne)
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
              
        diff=mass1-mass2
        
        T[i]=T2
        P[i]=P2/100000
        T_rat[i]=T2/T1
        P_rat[i]=P2/P1
        err[i]=diff
        
        P2=P2+stp
        
        
    
    fig,rT=plt.subplots() 
    rT.plot(T_rat, err)
    rT.set(xlabel='T2/T1', ylabel='Mass error (-)')
    rT.grid()
    
    fig,rT=plt.subplots() 
    rT.plot(P_rat, err)
    rT.set(xlabel='P2/P1', ylabel='Mass error (-)')
    rT.grid()
    
    return err, T_rat, P_rat
