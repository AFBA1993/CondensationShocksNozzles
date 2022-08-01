# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 20:26:27 2019

@author: USER
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from Coollib import Cool_PT, Cool_Psat, Cool_Tsat


def cndshock(P2,P1,T1,u1,q,f,EoS):
    prop1=Cool_PT(P1,T1,f,EoS)
    h1=prop1[3]
    rho1=prop1[1]
    propsat=Cool_Psat(P2,f,EoS)

    
#    psat=Cool_Tsat(T2,f,EoS)
#    psat=psat[0] 
    
    rho_g=1/propsat[2]
    h_g=propsat[3]
    

    h_l=propsat[7]
    
    X=(P1-P2+rho1*u1**2.0)*rho_g/(rho1*u1)**2
    y=1-X
      
    h2=X*h_g+y*h_l
    u2=(2*(h1-h2+0.5*u1**2.0+q))**0.5
    rho2=rho_g/X
    
          
    mass1=rho1*u1
    mass2=rho2*u2
          
    diff=mass1-mass2
    return diff



def cndshock_sol(P2,P1,T1,u1,q,f,EoS):
    prop1=Cool_PT(P1,T1,f,EoS)
    h1=prop1[3]
    rho1=prop1[1]
    propsat=Cool_Psat(P2,f,EoS)

    
#    psat=Cool_Tsat(T2,f,EoS)
#    psat=psat[0] 
    
    rho_g=1/propsat[2]
    h_g=propsat[3]
    

    h_l=propsat[7]
    
    X=(P1-P2+rho1*u1**2.0)*rho_g/(rho1*u1)**2
    y=1-X
      
    h2=X*h_g+y*h_l
    u2=(2*(h1-h2+0.5*u1**2.0+q))**0.5
    rho2=rho_g/X
              
    return rho2,X,u2





def cshock(P1,T1,u1,q,f,EoS):
    
    prop1=Cool_PT(P1,T1,f,EoS)
    h1=prop1[3]
    rho1=prop1[1]
    
    #Equations for cndshock
    P_trp=330*P1
    Pc=350*P1
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
        
        J_1[i]=(mass1)**2.0
        J_2[i]=-(P2-P1)/(1/rho2-1/rho1)
        
        h_1[i]=h1
        h_2[i]=(P2-P1)*(1/rho2+1/rho1)*0.5
        
        err[i]=diff 
        
        """
        
        if i>0:
            
            if  (err[i]*err[i-1])<0:
                if abs (err[i]-err[i-1])<3:
                    Pgss=P2-0.5*stp
                    c=c+1
                    if c==1:
                        P_wd=fsolve(cndshock,Pgss,(P1,T1,u1,q,f,EoS))
                        tol=cndshock(P_wd,P1,T1,u1,q,f,EoS)
                        if abs (tol)>1.8e+07:
                            print ("Pwd not found", tol)
                        wd_sol=cndshock_sol(P_wd,P1,T1,u1,q,f,EoS)
                        rho_wd=wd_sol[0]
                        X_wd=wd_sol[1]
                        u_wd=wd_sol[2]           
                        
                    
                        
                        
                        
                    if c==2:
                        P_sd=fsolve(cndshock,Pgss,(P1,T1,u1,q,f,EoS))
                        tol=cndshock(P_sd,P1,T1,u1,q,f,EoS)
                        if abs (tol)>1.8e+07:
                            print ("Psd not found", tol)
                        sd_sol=cndshock_sol(P_sd,P1,T1,u1,q,f,EoS)
                        rho_sd=wd_sol[0]
                        X_sd=sd_sol[1]
                        u_sd=sd_sol[2]   
                    
        """
        
        T[i]=T2
        P[i]=P2/100000
        T_rat[i]=T2/T1
        P_rat[i]=P2/P1
        
        
        
                    
        
        P2=P2+stp
        
        #print (P_wd)
      
        
        
    
    fig,rT=plt.subplots() 
    rT.plot(T_rat, err)
    rT.set(xlabel='T2/T1', ylabel='Mass error (-)')
    rT.grid()
    
    fig,rT=plt.subplots() 
    rT.plot(P_rat, err)
    rT.set(xlabel='P2/P1', ylabel='Mass error (-)')
    rT.grid()
    
    return err, T_rat, P_rat, J_1, J_2 #, P_wd, rho_wd, X_wd, u_wd, P_sd, rho_sd, X_sd, u_sd






