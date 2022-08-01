# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 19:13:33 2019

@author: USER
"""
import numpy as np
from scipy.optimize import fsolve
from Coollib import Cool_PT, Cool_Ps, Cool_hs 
from cflow_cool import  fnd_isenexp, fnd_throat,diff_isenexp
import matplotlib.pyplot as plt
from Cool_sat import p_sat, Ncl_rate
from scipy.integrate import quad
from scipy import interpolate

def Ix(x,P4,P3,A4,A3):
        m=(P4-P3)/(A4-A3)
        b=P3-m*A3
        y=m*x+b       
        return y

"""
def B2_Bier(x):
    if x<=2.01355789988531:
        y=0.008123539812761996*x**3.0+0.012445354905944056*x**2.0+0.009423166697206373*x+1.0000000000000002
    if x>2.01355789988531:
        y=-0.0010506010852403273*x**3.0+0.02646759857633423*x**2.0+0.025729978042863963*x+0.9852091434987051

    return y 
"""
def B2_Bier(x):
    x_data=np.zeros(4)
    y_data=np.zeros(4)
    
    x_data[0]=0    
    x_data[1]=2.25352643886249
    x_data[2]=4.28656558799125
    x_data[3]=6
    
    y_data[0]=1
    y_data[1]=1.16005380805109
    y_data[2]=1.4973850292927
    y_data[3]=1.86585703150537
    
    tck = interpolate.splrep(x_data, y_data)
    return interpolate.splev(x, tck)    

   
def BierB2_expansion(step,A2,P2,s2,f,EoS,Jcl,uplimT,xg):
    prop2=Cool_Ps(P2,s2,f,EoS,uplimT)
    T2=prop2[2]
    V2=1/prop2[1]
    h2=prop2[3]
    a2=prop2[6]
    u2=a2
    e2=h2+a2**2.0/2
    mass2=A2*a2/V2
    
    xth=0
    xi=0
    xf=6

    ne=10000
    L=xf-xi
    stp=L/ne
    Tc=304.1282
    
    
    Cool_rho=np.zeros(ne)
    Cool_x=np.zeros(ne)
    Cool_r=np.zeros(ne)
    Cool_P=np.zeros(ne)
    Cool_Pbar=np.zeros(ne)
    Cool_V=np.zeros(ne)
    Cool_T=np.zeros(ne)
    Cool_h=np.zeros(ne)
    Cool_s=np.zeros(ne)
    Cool_u=np.zeros(ne)
    Cool_M=np.zeros(ne)
    Cool_A=np.zeros(ne)
    Cool_S=np.zeros(ne)
    Cool_Jcl=np.zeros(ne)
    ei=np.zeros(ne)
    massi=np.zeros(ne)
    emass=np.zeros(ne)
    eener=np.zeros(ne)
    emom=np.zeros(ne)
    
    
    x=xi
    for i in range (ne):
        r=B2_Bier(x)
        Ai=A2*r
        
        if x==xth:            
            r=1
            Pi=P2
            si=s2
            propi=Cool_Ps(Pi,si,f,EoS,uplimT)
            ui=propi[6]
     
       
        if i>0:
          stpx=-0.001*Cool_P[i-1]
          Pi=fnd_isenexp(Cool_P[i-1],stpx,Ai,mass2,e2,s2,f,EoS,uplimT)
          
        si=s2  
        propi=Cool_Ps(Pi,si,f,EoS,uplimT)
        ui=(2*(h2+0.5*u2**2-propi[3]))**0.5                
               
        
        Mi=ui/propi[6]
        Cool_A[i]=Ai
        Cool_x[i]=x
        Cool_r[i]=r
        Cool_P[i]=propi[0]
        Cool_V[i]=1/propi[1]
        Cool_rho[i]=propi[1]
        Cool_T[i]=propi[2]
        #print (Cool_T[i])
        Cool_h[i]=propi[3]
        Cool_s[i]=propi[4]
        Cool_u[i]=ui
        Cool_M[i]=Mi
        Cool_Pbar[i]=Cool_P[i]/100000
        ei[i]=Cool_h[i]+0.5*Cool_u[i]**2.0
        massi[i]=Cool_rho[i]*Cool_u[i]*Cool_A[i] 
        
       
        
       
        if Cool_T[i]<Tc:
            Psat=p_sat(Cool_T[i],f,EoS)
            Cool_S[i]=Cool_P[i]/Psat
        
            if Cool_S[i]>1:
                Cool_Jcl[i]=Ncl_rate(Cool_P[i],Cool_T[i],f,EoS)
        if i>0: 
            
            I=I=quad(Ix,Cool_A[i-1],Cool_A[i],(Cool_P[i],Cool_P[i-1],Cool_A[i],Cool_A[i-1]))
            I=I[0]
            mom1=(Cool_P[i-1]+Cool_rho[i-1]*Cool_u[i-1]**2.0)*Cool_A[i-1]+I
            mom2=(Cool_P[i]+Cool_rho[i]*Cool_u[i]**2.0)*Cool_A[i]
            emom[i]=mom2-mom1    
            emass[i]=massi[i]-massi[i-1]
            eener[i]=ei[i]-ei[i-1]
        
        print(emom[i])
        
        
        
        x_cond=Cool_x[i-1]
        u_cond=Cool_u[i-1]
        T_cond=Cool_T[i-1]
        P_cond=Cool_P[i-1]
        J_cond=Cool_Jcl[i]
        rcrit_cond=0
        print(x,Cool_T[i])
        #if Cool_Jcl[i]>=Jcl:
            #break
        if x>=xg:
            break
        x+=stp
    """    
    plt.xlim(xi,x)
    plt.plot(Cool_x,Cool_P,Cool_S[i])
    """   
    return Cool_x,Cool_r,Cool_P,Cool_rho,Cool_T,Cool_h,Cool_s,Cool_u,Cool_M,Cool_S,Cool_Jcl,x_cond,u_cond,T_cond,P_cond,ne,Cool_A,J_cond,rcrit_cond, emom, emass, eener,ei






    

    
    
    
    
    



