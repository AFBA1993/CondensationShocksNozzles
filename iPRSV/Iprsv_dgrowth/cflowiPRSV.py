# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 15:46:46 2019

@author: sisea
"""
import numpy as np
from scipy.optimize import fsolve
from LibiPRSV import iPRSV_mas_PT, iPRSV_mas_Ps, iPRSV_mas_hs
import multiprocessing 

def pss_Shock (guess,step,k1,k2,k3,f):
  for i in range (1000):
    PT4=fsolve(diff_Shock,(guess,250),(k1,k2,k3,f))
    err=diff_Shock(PT4,k1,k2,k3,f)
    if abs(err[0])<=1.8e-08:
      if abs(err[1])<=1.8e-08:
        break
    guess=guess+step
  return PT4


def diff_thrt(P2,A1,A2,V1,h1,s1,f):
    s2=s1
    prop2=iPRSV_mas_Ps(f,P2,s2)
    V2=prop2[1]
    h2=prop2[3]
    a2=prop2[6]
    e2=h2+a2**2.0/2
    mass2=A2*a2/V2
    u1=(2*(e2-h1))**0.5
    mass1=A1*u1/V1
    y=mass1-mass2
    return y
  
def diff_isenexp(P3,A3,mass2,e2,s2,f):
  s3=s2
  Prop3=iPRSV_mas_Ps(f,P3,s3)
  V3=Prop3[1]
  h3=Prop3[3]
  u3=(2*(e2-h3))**0.5
  mass3=A3*u3/V3
  y=mass3-mass2
  return y
  
def diff_stagthr(P1,A1,A2,ho1,so1,f):
  s1=so1
  Propo1=iPRSV_mas_Ps(f,P1,s1)
  V1=Propo1[1]
  h1=Propo1[3]
  s1=Propo1[4]
  
  step=-0.005*P1
  guess=0.8*P1
  P2=fnd_throat(guess,step,A1,A2,V1,h1,s1,f)
  s2=s1
  
  Prop2=iPRSV_mas_Ps(f,P2,s2)
  
  a2=Prop2[6]
  h2=Prop2[3]
  e2=h2+a2**2.0/2
  u1=(2*(e2-h1))**0.5
  e1=h1+u1**2.0/2
  y=ho1-e1
  return y 

def diff_Shock(Z,k1,k2,k3,f):
  P4=Z[0]
  T4=Z[1]
  Prop4=iPRSV_mas_PT(f,P4,T4)
  V4=Prop4[1]
  h4=Prop4[3]
  u4=k1*V4
  P4f=k2-u4**2.0/V4
  h4f=k3-0.5*u4**2.0
  y=np.zeros(2)
  y[0]=P4-P4f
  y[1]=h4-h4f
  return y

def diff_Shockloc(r,P5r,A2,A5,P2,s2,mass2,e2,f,gss3,stp3,gss4,stp4,stp5):
  #Default guess value
  #gss3=P2, stp3=-0.01*P2; gss4=0.7*P1, stp4=-0.01*P1; gss5=P4, stp5=0.005*P2
  print(r)
  A3=r*A2
  P3=fnd_isenexp(gss3,stp3,A3,mass2,e2,s2,f)
  s3=s2
  Prop3=iPRSV_mas_Ps(f,P3,s3)
  V3=Prop3[1]
  h3=Prop3[3]
  a3=Prop3[6]
  u3=(2*(e2-h3))**0.5
  M3=u3/a3
  print('Upstream Normal Shock Wave achieved',M3)
  k1=u3/V3
  k2=P3+u3**2.0/V3
  k3=h3+u3**2.0/2
  PT4=pss_Shock(gss4,stp4,k1,k2,k3,f)  
  P4=PT4[0]
  T4=PT4[1]
  Prop4=iPRSV_mas_PT(f,P4,T4)
  V4=Prop4[1]
  h4=Prop4[3]
  s4=Prop4[4]
  a4=Prop4[6]
  u4=k1*V4
  M4=u4/a4
  e4=h4+u4**2.0/2
  mass4=u4*A3/V4
  print('Downstream Normal Shock Wave achieved',M4)
  P5=fnd_isenexp(P4,stp5,A5,mass4,e4,s4,f)
  s5=s4
  Prop5=iPRSV_mas_Ps(f,P5,s5)
  h5=Prop5[3]
  a5=Prop5[6]
  u5=(2*(e4-h5))**0.5
  M5=u5/a5
  print('Nozzle outlet achieved',M5)
  y=P5-P5r
  return y
 

def fnd_throat(guess,step,A1,A2,V1,h1,s1,f):
  error=np.zeros(1000)
  for i in range (1000):
    error[i]=diff_thrt(guess,A1,A2,V1,h1,s1,f)
    if i>=1:
      guess=guess+step
      error[i]=diff_thrt(guess,A1,A2,V1,h1,s1,f)
      if (error[i]*error[i-1])<0:
        guess=guess-0.5*step
        P2=fsolve(diff_thrt,guess,((A1,A2,V1,h1,s1,f)))
        err=diff_thrt(P2,A1,A2,V1,h1,s1,f)
        if abs(err)>1.8e-08:
          print ("Warning Throat pressure not found")
        break
  return P2


def fnd_stagthr(guessP1,A1,A2,ho1,so1,f):
  P1=fsolve(diff_stagthr,guessP1,(A1,A2,ho1,so1,f))
  error=diff_stagthr(P1,A1,A2,ho1,so1,f)
  if abs(error)>=1.8e-08:
    print ('from stagnation to nozzle error change the guess of the nozzle inlet')
  return P1

def fnd_isenexp(guess,step,A3,mass2,e2,s2,f):
  error=np.zeros(1000)
  for i in range (1000):
    error[i]=diff_isenexp(guess,A3,mass2,e2,s2,f)
    if i>=1:
      guess=guess+step
      error[i]=diff_isenexp(guess,A3,mass2,e2,s2,f)
      if (error[i]*error[i-1])<0:
        guess=guess-0.5*step
        P3=fsolve(diff_isenexp,guess,(A3,mass2,e2,s2,f))
        err=diff_isenexp(P3,A3,mass2,e2,s2,f)
        if abs(err)>1.8e-08:
            print ("Warning isentropic expansion pressure not found")
        break
  return P3 

def fnd_isenexp_try(inputs):
  guess=inputs[0]
  step=inputs[1]
  A3=inputs[2]
  mass2=inputs[3]
  e2=inputs[4]
  s2=inputs[5]
  f=inputs[6]
  error=np.zeros(1000)
  
  for i in range (1000):
    error[i]=diff_isenexp(guess,A3,mass2,e2,s2,f)
    if i>=1:
      guess=guess+step
      error[i]=diff_isenexp(guess,A3,mass2,e2,s2,f)
      if (error[i]*error[i-1])<0:
        guess=guess-0.5*step
        P3=fsolve(diff_isenexp,guess,(A3,mass2,e2,s2,f))
        err=diff_isenexp(P3,A3,mass2,e2,s2,f)
        if abs(err)>1.8e-08:
            print ("Warning isentropic expansion pressure not found")
        break
  return P3 















  
  
  
