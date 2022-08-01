# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 19:13:33 2019

@author: USER
"""
import numpy as np
from scipy.optimize import fsolve
from Coollib import Cool_PT, Cool_Ps, Cool_hs, Cool_Tsat, Cool_Psat
from cflow_cool import  fnd_isenexp, fnd_throat,pss_Shock, diff_Shockloc,diff_isenexp
import matplotlib.pyplot as plt
from Cool_sat import p_sat, Ncl_rate, T_sat
from scipy.integrate import quad
def Moore_nozzle(x):
  xth=2
  if x<=xth:
    Ar=-0.0625*x+1.127
  if x>xth:
    Ar=0.088*x+0.824
  return Ar


def Ix(x,P4,P3,A4,A3):
        m=(P4-P3)/(A4-A3)
        b=P3-m*A3
        y=m*x+b       
        return y



def Moore_expansion(step,A2,P2,s2,f,EoS,Jcl,xg,k):
  prop2=Cool_Ps(P2,s2,f,EoS)
  V2=1/prop2[1]
  h2=prop2[3]
  a2=prop2[6]
  u2=a2
  e2=h2+a2**2.0/2
  mass2=A2*a2/V2
  xth=2
  L=7
  div=62*10
  stp=1/div
  ne=div*L
     
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
  Cool_Sub=np.zeros(ne)
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
    if x<xth:
      Pi=fnd_isenexp(P2,0.07*P2,Ai,mass2,e2,s2,f,EoS)
      si=s2
      propi=Cool_Ps(Pi,si,f,EoS)
      ui=(2*(h2+0.5*u2**2-propi[3]))**0.5 
    if x==xth:
      r=1
      Pi=P2
      si=s2
      propi=Cool_Ps(Pi,si,f,EoS)
      ui=propi[6]
    if x>xth:
      Pi=fnd_isenexp(P2,step,Ai,mass2,e2,s2,f,EoS)
      si=s2
      propi=Cool_Ps(Pi,si,f,EoS)
      ui=(2*(h2+0.5*u2**2-propi[3]))**0.5
    
      
    if x>L:
      break
  
      
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
    
    if i>0:
        I=I=quad(Ix,Cool_A[i-1],Cool_A[i],(Cool_P[i],Cool_P[i-1],Cool_A[i],Cool_A[i-1]))
        I=I[0]
        mom1=(Cool_P[i-1]+Cool_rho[i-1]*Cool_u[i-1]**2.0)*Cool_A[i-1]+I
        mom2=(Cool_P[i]+Cool_rho[i]*Cool_u[i]**2.0)*Cool_A[i]
        emom[i]=mom2-mom1    
        emass[i]=massi[i]-massi[i-1]
        eener[i]=ei[i]-ei[i-1]
    
    Cool_Sub[i]=T_sat(Cool_P[i],f,EoS)-Cool_T[i]
    
    
    Cool_S[i]=Cool_P[i]/p_sat(Cool_T[i],f,EoS)
    if Cool_S[i]>1:
        Ji=Ncl_rate(Cool_P[i],Cool_T[i],f,EoS,k)
        Cool_Jcl[i]=Ji[0]
        rcrit=Ji[1]
        Cool_rcrit=Ji[1]
    
    
    print(Cool_x[i],Cool_S[i],Cool_Jcl[i])
    
    x_cond=Cool_x[i]
    u_cond=Cool_u[i]
    T_cond=Cool_T[i]
    P_cond=Cool_P[i]
    J_cond=Cool_Jcl[i]
    rcrit_cond=0
    
    if f=="Water":
      if Cool_T[i]<=273.16:       
         break
     
    if Cool_x[i]>=xg:
        break
    if Cool_Jcl[i]>0:
        break
         
   
    x=x+stp

  return Cool_x,Cool_r,Cool_P,Cool_rho,Cool_T,Cool_h,Cool_s,Cool_u,Cool_M,Cool_S,Cool_Jcl,x_cond,u_cond,T_cond,P_cond,ne,Cool_A,J_cond,rcrit_cond,Cool_Sub,Cool_rcrit,emom, emass, eener,ei




