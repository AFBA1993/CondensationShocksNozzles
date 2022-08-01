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

def Moore_nozzle(x):
  xth=2
  if x<=xth:
    Ar=-0.0625*x+1.127
  if x>=xth:
    Ar=0.088*x+0.824
  return Ar

def Moore_expansion(step,A2,P2,s2,f):
  prop2=iPRSV_mas_Ps(f,P2,s2)
  V2=prop2[1]
  h2=prop2[3]
  a2=prop2[6]
  u2=a2
  e2=h2+a2**2.0/2
  mass2=A2*a2/V2
  xth=2
  L=7
  div=64
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
   
  x=0
   
  for i in range (ne):
    r=Moore_nozzle(x)
    Ai=r*A2
    if x<xth:
      Pi=fnd_isenexp(P2,0.07*P2,Ai,mass2,e2,s2,f)
      si=s2
      propi=iPRSV_mas_Ps(f,Pi,si)
      ui=(2*(h2+0.5*u2**2-propi[3]))**0.5 
    if x==xth:
      r=1
      Pi=P2
      si=s2
      propi=iPRSV_mas_Ps(f,Pi,si)
      ui=propi[6]
    if x>xth:
      Pi=fnd_isenexp(P2,step,Ai,mass2,e2,s2,f)
      si=s2
      propi=iPRSV_mas_Ps(f,Pi,si)
      ui=(2*(h2+0.5*u2**2-propi[3]))**0.5
    if x>L:
      break
      
    Mi=ui/propi[6]
    iPRSV_A[i]=Ai
    iPRSV_x[i]=x
    iPRSV_r[i]=r
    iPRSV_P[i]=propi[0]/100000
    iPRSV_V[i]=propi[1]
    iPRSV_rho[i]=1/propi[1]
    iPRSV_T[i]=propi[2]
    iPRSV_h[i]=propi[3]
    iPRSV_s[i]=propi[4]
    iPRSV_u[i]=ui
    iPRSV_M[i]=Mi
    print(Mi)
    
    if f=="Water":
      if iPRSV_T[i]<=273.16:
       if iPRSV_P[i]<=(611.6548008968684/100000):
         break  
    
    x=x+stp
    
  ig,rT=plt.subplots() 
  rT.plot(iPRSV_x, iPRSV_T)
  rT.set(xlabel='Lenght (-)', ylabel='Temperature (K)')
  rT.grid() 
  return iPRSV_x,iPRSV_r,iPRSV_P,iPRSV_rho,iPRSV_T,iPRSV_h,iPRSV_s,iPRSV_u,iPRSV_M





