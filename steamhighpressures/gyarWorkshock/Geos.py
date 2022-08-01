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
from Cool_sat import p_sat, Ncl_rate, Ncl_rate2
from scipy.integrate import quad

def Ix(x,P4,P3,A4,A3):
        m=(P4-P3)/(A4-A3)
        b=P3-m*A3
        y=m*x+b       
        return y

def B4M(x):
    coeff=np.zeros(12)
    coeff[0]=1#0.99998
    coeff[1]=-0.000364615
    coeff[2]=0.00642422
    coeff[3]=-0.000672132
    coeff[4]=3.36339e-05
    coeff[5]=9.4511e-07
    coeff[6]=-2.03138e-07
    coeff[7]=1.01671e-08
    coeff[8]=-2.2708e-10
    coeff[9]=1.95149e-12
    coeff[10]=0
    coeff[11]=0
    y=coeff[0]+coeff[1]*x+coeff[2]*x**2.0+coeff[3]*x**3.0+coeff[4]*x**4+coeff[5]*x**5.0+coeff[6]*x**6.0+coeff[7]*x**7.0+coeff[8]*x**8.0+coeff[9]*x**9.0+coeff[10]*x**10+coeff[11]*x**11
    return y

def B2M(x):
    coeff=np.zeros(12)
    coeff[0]=1#0.9999942992
    coeff[1]=1.922629195e-05
    coeff[2]=0.0002620341402
    coeff[3]=-6.198438881e-06
    coeff[4]=4.392952345e-09
    coeff[5]=6.293312484e-09
    coeff[6]=-2.176828098e-10
    coeff[7]=3.402165174e-12
    coeff[8]=-2.427327457e-14
    coeff[9]=1.939688604e-17
    coeff[10]=6.984506287e-19
    coeff[11]=-2.804815032e-21    
    y=coeff[0]+coeff[1]*x+coeff[2]*x**2.0+coeff[3]*x**3.0+coeff[4]*x**4+coeff[5]*x**5.0+coeff[6]*x**6.0+coeff[7]*x**7.0+coeff[8]*x**8.0+coeff[9]*x**9.0+coeff[10]*x**10+coeff[11]*x**11
    return y



def B2M_expansion(step,A2,P2,s2,f,EoS,Jcl,To,lim):
    prop2=Cool_Ps(P2,s2,f,EoS,To,lim)
    V2=1/prop2[1]
    h2=prop2[3]
    a2=prop2[6]
    u2=a2
    e2=h2+a2**2.0/2
    mass2=A2*a2/V2
    
    gssT=To
    
    xth=0
    xi=0
    xf=100
    ne=800
    L=xf-xi
    
    stp=L/ne

    Cool_rho=np.zeros(ne)
    Cool_x=np.zeros(ne)
    Cool_r=np.zeros(ne)
    Cool_P=np.zeros(ne)
    Cool_V=np.zeros(ne)
    Cool_T=np.zeros(ne)
    Cool_h=np.zeros(ne)
    Cool_s=np.zeros(ne)
    Cool_u=np.zeros(ne)
    Cool_M=np.zeros(ne)
    Cool_A=np.zeros(ne)
    Cool_S=np.zeros(ne)
    Cool_Jcl=np.zeros(ne)
    
    x=xi
    
    for i in range (ne):
        r=B2M(x)
        Ai=A2*r
        if i==0:
            gssT=To
       
        
        if x==xth:
            gssT=To
            r=1
            Pi=P2
            si=s2
            propi=Cool_Ps(Pi,si,f,EoS,gssT,lim)
            ui=propi[6]
        if x>xth:
          lim=Cool_T[i-1]-5
          Pi=fnd_isenexp(Cool_P[i-1],-0.0001*Cool_P[i-1],Ai,mass2,e2,s2,f,EoS,Cool_T[i-1],lim)
          si=s2
          propi=Cool_Ps(Pi,si,f,EoS,Cool_T[i-1],lim)
          ui=(2*(h2+0.5*u2**2-propi[3]))**0.5
          
        Mi=ui/propi[6]
        Cool_A[i]=Ai
        Cool_x[i]=x
        Cool_r[i]=r
        Cool_P[i]=propi[0]
        Cool_V[i]=1/propi[1]
        Cool_rho[i]=propi[1]
        Cool_T[i]=propi[2]
        Cool_h[i]=propi[3]
        Cool_s[i]=propi[4]
        Cool_u[i]=ui
        Cool_M[i]=Mi
        
        Cool_S[i]=Cool_P[i]/p_sat(Cool_T[i],f,EoS)
        if Cool_S[i]>1:
          Cool_Jcl[i]=Ncl_rate2(Cool_P[i],Cool_T[i],f,EoS)
        
        
        
        x_cond=Cool_x[i]
        u_cond=Cool_u[i]
        T_cond=Cool_T[i]
        P_cond=Cool_P[i]
        
        if f=="Water":
          if Cool_T[i]<=273.16:       
             break
        if Cool_Jcl[i]>=Jcl:
            break
        if x>=60:
           break
            
        print (x, Cool_T[i], Cool_Jcl[i])
        x=x+stp
       
        
    fig,rT=plt.subplots() 
    rT.plot(Cool_x, Cool_T)
    rT.set(xlabel='Lenght (-)', ylabel='Temperature (K)')
    rT.grid() 
    return Cool_x,Cool_r,Cool_P,Cool_rho,Cool_T,Cool_h,Cool_s,Cool_u,Cool_M,Cool_S,Cool_Jcl,x_cond,u_cond,T_cond,P_cond,ne,Cool_A

def B4M_expansion(step,A2,P2,s2,f,EoS,Jcl,To,lim,xg):
    prop2=Cool_Ps(P2,s2,f,EoS,To,lim)
    V2=1/prop2[1]
    h2=prop2[3]
    a2=prop2[6]
    u2=a2
    e2=h2+a2**2.0/2
    mass2=A2*a2/V2
    
    gssT=To
    
    xth=0
    xi=0
    xf=30
    ne=10000
    L=xf-xi
    
    stp=L/ne

    Cool_rho=np.zeros(ne)
    Cool_x=np.zeros(ne)
    Cool_r=np.zeros(ne)
    Cool_P=np.zeros(ne)
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
    gssP=P2
    stpg=-0.001*P2
    
    jj=0
    for i in range (ne):
        r=B4M(x)
        Ai=A2*r
        if i==0:
            gssT=To
            
       
        
        if x==xth:
            gssT=To
            r=1
            Pi=P2
            si=s2
            propi=Cool_Ps(Pi,si,f,EoS,gssT,lim)
            ui=propi[6]
        if x>xth:
          #gsssT=Cool_T[0]
          #lim=Cool_T[i-1]-30
          if x>0:
              #if jj==0:
                #x=0.5
                #r=B4M(x)
                #Ai=A2*r
                #jj=1
              lim=Cool_T[i-1]-5
              gssP=1*Cool_P[i-1]
              stpg=-0.001*Cool_P[i-1]
              lim=Cool_T[i-1]-5
              gsssT=Cool_T[i-1]
              #gssP=Cool_P[0]
              #stpg=-0.001*Cool_P[0]
              #gsssT=Cool_T[0]+30
              #lim=Cool_T[i-1]-5
          Pi=fnd_isenexp(gssP,stpg,Ai,mass2,e2,s2,f,EoS,gsssT,lim)
          si=s2
          propi=Cool_Ps(Pi,si,f,EoS,Cool_T[i-1],lim)
          ui=(2*(h2+0.5*u2**2-propi[3]))**0.5
          
        Mi=ui/propi[6]
        Cool_A[i]=Ai
        Cool_x[i]=x
        Cool_r[i]=r
        Cool_P[i]=propi[0]
        Cool_V[i]=1/propi[1]
        Cool_rho[i]=propi[1]
        Cool_T[i]=propi[2]
        Cool_h[i]=propi[3]
        Cool_s[i]=propi[4]
        Cool_u[i]=ui
        Cool_M[i]=Mi
        ei[i]=Cool_h[i]+0.5*Cool_u[i]**2.0
        massi[i]=Cool_rho[i]*Cool_u[i]*Cool_A[i] 
        
        Cool_S[i]=Cool_P[i]/p_sat(Cool_T[i],f,EoS)
        if Cool_S[i]>1:
          Cool_Jcl[i]=Ncl_rate2(Cool_P[i],Cool_T[i],f,EoS)
        
        if i>0: 
            
            I=I=quad(Ix,Cool_A[i-1],Cool_A[i],(Cool_P[i],Cool_P[i-1],Cool_A[i],Cool_A[i-1]))
            I=I[0]
            mom1=(Cool_P[i-1]+Cool_rho[i-1]*Cool_u[i-1]**2.0)*Cool_A[i-1]+I
            mom2=(Cool_P[i]+Cool_rho[i]*Cool_u[i]**2.0)*Cool_A[i]
            emom[i]=mom2-mom1    
            emass[i]=massi[i]-massi[i-1]
            eener[i]=ei[i]-ei[i-1]
        
        
        
        x_cond=Cool_x[i]
        u_cond=Cool_u[i]
        T_cond=Cool_T[i]
        P_cond=Cool_P[i]
        J_cond=Cool_Jcl[i]
        rcrit_cond=0
        
        if f=="Water":
          if Cool_T[i]<=273.16:       
             break
        #if Cool_Jcl[i]>=Jcl:
            #break
        if x>=xg:
           break
        if x>=2.0:
            stp=L/ne
        else:
            stp=L/400
            emom[i]=1e-08
        print (x, emom[i])
        x=x+stp
       
      
    fig,rT=plt.subplots() 
    rT.plot(Cool_x, Cool_P)
    rT.set(xlabel='Lenght (-)', ylabel='Temperature (K)')
    rT.grid()
    
    return Cool_x,Cool_r,Cool_P,Cool_rho,Cool_T,Cool_h,Cool_s,Cool_u,Cool_M,Cool_S,Cool_Jcl,x_cond,u_cond,T_cond,P_cond,ne,Cool_A,J_cond,rcrit_cond, emom, emass, eener,ei





