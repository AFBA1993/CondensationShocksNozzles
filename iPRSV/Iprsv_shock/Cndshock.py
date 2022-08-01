# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 09:50:46 2019

@author: sisea
"""
from LibiPRSV import iPRSV_mas_Ps, iPRSV_mas_PT, iPRSV_mas_PT_sat,iPRSV_mas_VT
from LibiPRSVsat import fnd_psatiPRSV, fnd_Tsat
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt





def cndshock(P2,P1,T1,u1,q,f):
    prop1=iPRSV_mas_PT(f,P1,T1)
    h1=prop1[3]
    rho1=1/prop1[1]
    
    T2=fnd_Tsat(P2,f)
    propsat2=iPRSV_mas_PT_sat(f,P2,T2)
    
    rho_g=1/propsat2[2]
    h_g=propsat2[3]
    h_l=propsat2[7]
        
    X=(P1-P2+rho1*u1**2.0)*rho_g/(rho1*u1)**2
    y=1-X
      
    h2=X*h_g+y*h_l
    u2=(2*(h1-h2+0.5*u1**2.0+q))**0.5
    rho2=rho_g/X
    
      
    mass1=rho1*u1
    mass2=rho2*u2
      
    diff=mass1-mass2
    
    return diff



def cndshock_sol(P2,P1,T1,u1,q,f):
    prop1=iPRSV_mas_PT(f,P1,T1)
    h1=prop1[3]
    rho1=1/prop1[1]
    
    T2=fnd_Tsat(P2,f)
    propsat2=iPRSV_mas_PT_sat(f,P2,T2)
    
    rho_g=1/propsat2[2]
    h_g=propsat2[3]
    h_l=propsat2[7]
        
    X=(P1-P2+rho1*u1**2.0)*rho_g/(rho1*u1)**2
    y=1-X
      
    h2=X*h_g+y*h_l
    u2=(2*(h1-h2+0.5*u1**2.0+q))**0.5
    rho2=rho_g/X
    
    return rho2,X,u2
    

def cshock(P1,T1,u1,q,f):
    prop1=iPRSV_mas_PT(f,P1,T1)
    h1=prop1[3]
    rho1=1/prop1[1]
    P_trp=1*P1
    Pc=2.2*P1
    Pdiff=Pc-P_trp
    ne=1000
    
    err=np.zeros(ne)
    T=np.zeros(ne)
    P=np.zeros(ne)
    T_rat=np.zeros(ne)
    P_rat=np.zeros(ne)
    J_1=np.zeros(ne)
    J_2=np.zeros(ne)
    h_1=np.zeros(ne)
    h_2=np.zeros(ne)
    
    stp=Pdiff/ne
    P2=P_trp
    c=0
    
    for i in range (ne):
        
        T2=fnd_Tsat(P2,f)
        propsat2=iPRSV_mas_PT_sat(f,P2,T2)
        
        rho_g=1/propsat2[2]
        h_g=propsat2[3]
        h_l=propsat2[7]
            
        X=(P1-P2+rho1*u1**2.0)*rho_g/(rho1*u1)**2
        y=1-X
          
        h2=X*h_g+y*h_l
        u2=(2*(h1-h2+0.5*u1**2.0+q))**0.5
        rho2=rho_g/X
        
          
        mass1=rho1*u1
        mass2=rho2*u2
          
        diff=mass1-mass2
        
        J_1[i]=(mass1)**2.0
        J_2[i]=-(P2-P1)/(1/rho2-1/rho1)
        
        h_1[i]=h1
        h_2[i]=(P2-P1)*(1/rho2+1/rho1)*0.5
        
        err[i]=diff
        
        
        
        if i>0:
                
                if  (err[i]*err[i-1])<0:
                    if abs (err[i]-err[i-1])<3:
                        Pgss=P2-0.5*stp
                        c=c+1
                        if c==1:
                            P_wd=fsolve(cndshock,Pgss,(P1,T1,u1,q,f))
                            tol=cndshock(P_wd,P1,T1,u1,q,f)
                            if abs (tol)>1.8e+07:
                                print ("Pwd not found", tol)
                            wd_sol=cndshock_sol(P_wd,P1,T1,u1,q,f)
                            rho_wd=wd_sol[0]
                            X_wd=wd_sol[1]
                            u_wd=wd_sol[2]
                            
                        if c==2:
                            P_sd=fsolve(cndshock,Pgss,(P1,T1,u1,q,f))
                            tol=cndshock(P_sd,P1,T1,u1,q,f)
                            if abs (tol)>1.8e+07:
                                print ("Psd not found", tol)
                            sd_sol=cndshock_sol(P_sd,P1,T1,u1,q,f)
                            rho_sd=sd_sol[0]
                            X_sd=sd_sol[1]
                            u_sd=sd_sol[2]  
        
        
        T[i]=T2
        P[i]=P2/100000
        T_rat[i]=T2/T1
        P_rat[i]=P2/P1
        
        
        P2=P2+stp
        
    fig,rT=plt.subplots() 
    rT.plot(T_rat, err)
    rT.set(xlabel='T2/T1', ylabel='Mass error (-)')
    rT.grid()
    
    fig,rT=plt.subplots() 
    rT.plot(P_rat, err)
    rT.set(xlabel='P2/P1', ylabel='Mass error (-)')
    rT.grid()
        
        
        
    return err, T_rat, P_rat, J_1, J_2, P_wd, rho_wd, X_wd, u_wd, P_sd, rho_sd, X_sd, u_sd
    
    
    
"""
f="Water"    
P3=9454.815732892805
T3=281.8525302855678
u3=534.1874275918426

cnd=cshock(P3,T3,u3,0,f)
"""














