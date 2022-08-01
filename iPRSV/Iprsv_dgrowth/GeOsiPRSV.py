# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 18:24:41 2019

@author: sisea
"""
import numpy as np
from scipy.optimize import fsolve
import math
from LibiPRSV import iPRSV_mas_Ps
from cflowiPRSV import fnd_isenexp
import matplotlib.pyplot as plt
from LibiPRSVsat import fnd_psatiPRSV, Ncl_rate, fnd_Tsat,fnd_psatiPRSV_thread,fnd_psatiPRSV_Process
from scipy.integrate import quad
from scipy import interpolate

def Ix(x,P4,P3,A4,A3):
        m=(P4-P3)/(A4-A3)
        b=P3-m*A3
        y=m*x+b       
        return y

def Moore_nozzle(x):
  xth=2
  if x<=xth:
    Ar=-0.0625*x+1.127
  if x>xth:
    Ar=0.088*x+0.824
  return Ar


def f1(x):
    x_points = [ 2, 2.3399999999999998579, 2.8020000000000000462, 3.3000000000000002665, 3.7400000000000002132, 4.3419999999999996376, 5.3390000000000004121]
    y_points = [4,3.404,2.75,2.188,1.784,1.406,1.086]

    tck = interpolate.splrep(x_points, y_points,k=2)
    return interpolate.splev(x, tck)
def f2(x):
    x_points = [ 5.3390000000000004121,8.22]
    y_points = [1.086, 1]

    tck = interpolate.splrep(x_points, y_points,k=1)
    return interpolate.splev(x, tck)
def f3(x):
    x_points = [8.22, 9.34, 10.188, 11.104, 11.96, 13.398,16]
    y_points = [1,1.024, 1.064, 1.128, 1.22,1.414,1.94186]

    tck = interpolate.splrep(x_points, y_points)
    return interpolate.splev(x, tck)

def Moses_nozzle(x):
    if (x>=2) & (x<=5.3390000000000004121):
        y=f1(x)
    if (x>5.3390000000000004121) & (x<=8.22):
        y=f2(x)
    if (x>8.22):
        y=f3(x)
    return y



def Moore_expansion(step,grid,A2,P2,s2,f,k,xg):
  prop2=iPRSV_mas_Ps(f,P2,s2)
  V2=prop2[1]
  h2=prop2[3]
  a2=prop2[6]
  u2=a2
  e2=h2+a2**2.0/2
  mass2=A2*a2/V2
  xth=2
  L=7
  div=grid
  stp=1/div
  ne=div*L
     
  iPRSV_rho=np.zeros(ne)
  iPRSV_x=np.zeros(ne)
  iPRSV_r=np.zeros(ne)
  iPRSV_P=np.zeros(ne)
  iPRSV_V=np.zeros(ne)
  iPRSV_T=np.zeros(ne)
  iPRSV_h=np.zeros(ne)
  iPRSV_s=np.zeros(ne)
  iPRSV_u=np.zeros(ne)
  iPRSV_M=np.zeros(ne)
  iPRSV_A=np.zeros(ne)
  iPRSV_S=np.zeros(ne)
  iPRSV_Subc=np.zeros(ne)
  iPRSV_psat=np.zeros(ne)
  iPRSV_Tsat=np.zeros(ne)
  iPRSV_S=np.zeros(ne)
  iPRSV_Jcl=np.zeros(ne)
  iPRSV_Z=np.zeros(ne)
  Cool_rcrit=np.zeros(ne)
  ei=np.zeros(ne)
  massi=np.zeros(ne)
  emass=np.zeros(ne)
  eener=np.zeros(ne)
  emom=np.zeros(ne)
  
  
  
  x=0
   
  for i in range (ne):
    r=Moore_nozzle(x)
    Ai=r*A2
    if i==0:
        if x<xth:
          Pi=fnd_isenexp(P2,0.07*P2,Ai,mass2,e2,s2,f)
          si=s2
          propi=iPRSV_mas_Ps(f,Pi,si)
          ui=(2*(h2+0.5*u2**2-propi[3]))**0.5 
    if i>0:
        Pi=fnd_isenexp(iPRSV_P[i-1],-0.005*iPRSV_P[i-1],Ai,mass2,e2,s2,f)
        si=s2
        propi=iPRSV_mas_Ps(f,Pi,si)
        ui=(2*(h2+0.5*u2**2-propi[3]))**0.5 
          
    if x==xth:
      r=1
      Pi=P2
      si=s2
      propi=iPRSV_mas_Ps(f,Pi,si)
      ui=propi[6]
    """
    if x>xth:
      Pi=fnd_isenexp(P2,step,Ai,mass2,e2,s2,f)
      si=s2
      propi=iPRSV_mas_Ps(f,Pi,si)
      ui=(2*(h2+0.5*u2**2-propi[3]))**0.5
    """
    if x>L:
      break
      
    Mi=ui/propi[6]
    iPRSV_A[i]=Ai
    iPRSV_x[i]=x
    iPRSV_r[i]=r
    iPRSV_P[i]=propi[0]
    iPRSV_V[i]=propi[1]
    iPRSV_rho[i]=1/propi[1]
    iPRSV_T[i]=propi[2]
    iPRSV_h[i]=propi[3]
    iPRSV_s[i]=propi[4]
    iPRSV_Z[i]=propi[5]
    iPRSV_u[i]=ui
    iPRSV_M[i]=Mi

    
    iPRSV_psat[i]=fnd_psatiPRSV_Process(f,iPRSV_T[i])
    #iPRSV_psat[i]=fnd_psatiPRSV(f,iPRSV_T[i])
    
    iPRSV_Tsat[i]=fnd_Tsat(iPRSV_P[i],f)
    iPRSV_Subc[i]=iPRSV_Tsat[i]-iPRSV_T[i]
    
    iPRSV_S[i]=iPRSV_P[i]/iPRSV_psat[i]
    Jcl=Ncl_rate(f,iPRSV_P[i],iPRSV_T[i],iPRSV_psat[i],iPRSV_Tsat[i],k)
    iPRSV_Jcl[i]=Jcl[0]
    Cool_rcrit[i]=Jcl[1]
    if i>0:
        I=quad(Ix,iPRSV_A[i-1],iPRSV_A[i],(iPRSV_P[i],iPRSV_P[i-1],iPRSV_A[i],iPRSV_A[i-1]))
        I=I[0]
        mom1=(iPRSV_P[i-1]+iPRSV_rho[i-1]*iPRSV_u[i-1]**2.0)*iPRSV_A[i-1]+I
        mom2=(iPRSV_P[i]+iPRSV_rho[i]*iPRSV_u[i]**2.0)*iPRSV_A[i]
        emom[i]=mom2-mom1
        emass[i]=massi[i]-massi[i-1]
        eener[i]=ei[i]-ei[i-1]
        print(x,"errormom", emom[i])
    
    
    if f=="Water":
      if iPRSV_T[i]<=273.16:
         break 
      if iPRSV_P[i]<=(611.6548008968684/100000):
         break  
    
    
    
    x_cond=iPRSV_x[i]
    u_cond=iPRSV_u[i]
    T_cond=iPRSV_T[i]
    P_cond=iPRSV_P[i]
    Jcl_cond=iPRSV_Jcl[i]
    rcrit_cond=0
    
    if iPRSV_Jcl[i]>0:
        break
    #print(x,iPRSV_S[i],iPRSV_Jcl[i])

    x=x+stp
  """ 
  ig,rT=plt.subplots() 
  rT.plot(iPRSV_x, iPRSV_T)
  rT.set(xlabel='Lenght (-)', ylabel='Temperature (K)')
  rT.grid()
  """
  
  return iPRSV_x,iPRSV_r,iPRSV_P,iPRSV_rho,iPRSV_T,iPRSV_h,iPRSV_s,iPRSV_u,iPRSV_M,iPRSV_S,iPRSV_Jcl,x_cond,u_cond,T_cond,P_cond,ne,iPRSV_A, Jcl_cond,rcrit_cond,iPRSV_Subc,Cool_rcrit,emom,emass,eener,ei  


def Moses_expansion(step,grid,A2,P2,s2,f,k,xg):
  prop2=iPRSV_mas_Ps(f,P2,s2)
  V2=prop2[1]
  h2=prop2[3]
  a2=prop2[6]
  u2=a2
  e2=h2+a2**2.0/2
  mass2=A2*a2/V2
  xth=8.22
  L=16-6
  div=grid
  stp=L/div
  ne=grid
     
  iPRSV_rho=np.zeros(ne)
  iPRSV_x=np.zeros(ne)
  iPRSV_r=np.zeros(ne)
  iPRSV_P=np.zeros(ne)
  iPRSV_V=np.zeros(ne)
  iPRSV_T=np.zeros(ne)
  iPRSV_h=np.zeros(ne)
  iPRSV_s=np.zeros(ne)
  iPRSV_u=np.zeros(ne)
  iPRSV_M=np.zeros(ne)
  iPRSV_A=np.zeros(ne)
  iPRSV_S=np.zeros(ne)
  iPRSV_Subc=np.zeros(ne)
  iPRSV_psat=np.zeros(ne)
  iPRSV_Tsat=np.zeros(ne)
  iPRSV_S=np.zeros(ne)
  iPRSV_Jcl=np.zeros(ne)
  iPRSV_Z=np.zeros(ne)
  Cool_rcrit=np.zeros(ne)
  ei=np.zeros(ne)
  massi=np.zeros(ne)
  emass=np.zeros(ne)
  eener=np.zeros(ne)
  emom=np.zeros(ne)
  
  
  
  x=6
   
  for i in range (ne):
    r=Moses_nozzle(x)
    Ai=r*A2
    if i==0:
        if x<xth:
          Pi=fnd_isenexp(P2,0.07*P2,Ai,mass2,e2,s2,f)
          si=s2
          propi=iPRSV_mas_Ps(f,Pi,si)
          ui=(2*(h2+0.5*u2**2-propi[3]))**0.5 
    if i>0:
        Pi=fnd_isenexp(iPRSV_P[i-1],-0.0001*iPRSV_P[i-1],Ai,mass2,e2,s2,f)
        si=s2
        propi=iPRSV_mas_Ps(f,Pi,si)
        ui=(2*(h2+0.5*u2**2-propi[3]))**0.5 
          
    if x==xth:
      r=1
      Pi=P2
      si=s2
      propi=iPRSV_mas_Ps(f,Pi,si)
      ui=propi[6]
    """
    if x>xth:
      Pi=fnd_isenexp(P2,step,Ai,mass2,e2,s2,f)
      si=s2
      propi=iPRSV_mas_Ps(f,Pi,si)
      ui=(2*(h2+0.5*u2**2-propi[3]))**0.5
    """
    #if x>L:
      #break
    print(i)
    Mi=ui/propi[6]
    iPRSV_A[i]=Ai
    iPRSV_x[i]=x
    iPRSV_r[i]=r
    iPRSV_P[i]=propi[0]
    iPRSV_V[i]=propi[1]
    iPRSV_rho[i]=1/propi[1]
    iPRSV_T[i]=propi[2]
    iPRSV_h[i]=propi[3]
    iPRSV_s[i]=propi[4]
    iPRSV_Z[i]=propi[5]
    iPRSV_u[i]=ui
    iPRSV_M[i]=Mi

    
    iPRSV_psat[i]=fnd_psatiPRSV_Process(f,iPRSV_T[i])
    #iPRSV_psat[i]=fnd_psatiPRSV(f,iPRSV_T[i])
    
    iPRSV_Tsat[i]=fnd_Tsat(iPRSV_P[i],f)
    iPRSV_Subc[i]=iPRSV_Tsat[i]-iPRSV_T[i]
    
    iPRSV_S[i]=iPRSV_P[i]/iPRSV_psat[i]
    Jcl=Ncl_rate(f,iPRSV_P[i],iPRSV_T[i],iPRSV_psat[i],iPRSV_Tsat[i],k)
    iPRSV_Jcl[i]=Jcl[0]
    Cool_rcrit[i]=Jcl[1]
    if i>0:
        I=quad(Ix,iPRSV_A[i-1],iPRSV_A[i],(iPRSV_P[i],iPRSV_P[i-1],iPRSV_A[i],iPRSV_A[i-1]))
        I=I[0]
        mom1=(iPRSV_P[i-1]+iPRSV_rho[i-1]*iPRSV_u[i-1]**2.0)*iPRSV_A[i-1]+I
        mom2=(iPRSV_P[i]+iPRSV_rho[i]*iPRSV_u[i]**2.0)*iPRSV_A[i]
        emom[i]=mom2-mom1
        emass[i]=massi[i]-massi[i-1]
        eener[i]=ei[i]-ei[i-1]
        print("errormom", emom[i])
    
    
    if f=="Water":
      if iPRSV_T[i]<=273.16:
         break 
      if iPRSV_P[i]<=(611.6548008968684/100000):
         break  
    
    
    
    x_cond=iPRSV_x[i]
    u_cond=iPRSV_u[i]
    T_cond=iPRSV_T[i]
    P_cond=iPRSV_P[i]
    Jcl_cond=iPRSV_Jcl[i]
    rcrit_cond=0
    
    if x>=xg:
        break
    
    if iPRSV_Jcl[i]>0:
        break
    print(x,emom[i])

    x=x+stp
    
  return iPRSV_x,iPRSV_r,iPRSV_P,iPRSV_rho,iPRSV_T,iPRSV_h,iPRSV_s,iPRSV_u,iPRSV_M,iPRSV_S,iPRSV_Jcl,x_cond,u_cond,T_cond,P_cond,ne,iPRSV_A, Jcl_cond,rcrit_cond,iPRSV_Subc,Cool_rcrit,emom,emass,eener,ei 






    