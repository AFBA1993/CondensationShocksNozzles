# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 22:06:51 2019

@author: USER
"""
import numpy as np
from scipy.optimize import fsolve
from Coollib import Cool_PT, Cool_Ps

def diff_thrt(P2,A1,A2,V1,h1,s1,f,EoS):
    s2=s1
    prop2=Cool_Ps(P2,s2,f,EoS)
    V2=1/prop2[1]
    h2=prop2[3]
    a2=prop2[6]
    e2=h2+a2**2.0/2
    mass2=A2*a2/V2
    u1=(2*(e2-h1))**0.5
    mass1=A1*u1/V1
    y=mass1-mass2
    return y

def diff_stagthr(P1,A1,A2,ho1,so1,f,EoS):
  s1=so1
  Propo1=Cool_Ps(P1,s1,f,EoS)
  V1=1/Propo1[1]
  h1=Propo1[3]
  s1=Propo1[4]
  
  step=-0.01*P1
  guess=0.8*P1
  P2=fnd_throat(guess,step,A1,A2,V1,h1,s1,f,EoS)
  s2=s1
  
  Prop2=Cool_Ps(P2,s2,f,EoS)
  
  a2=Prop2[6]
  h2=Prop2[3]
  e2=h2+a2**2.0/2
  u1=(2*(e2-h1))**0.5
  e1=h1+u1**2.0/2
  y=ho1-e1
  return y


def fnd_throat(guess,step,A1,A2,V1,h1,s1,f,EoS):
  error=np.zeros(1000)
  for i in range (1000):
    print(guess)
    error[i]=diff_thrt(guess,A1,A2,V1,h1,s1,f,EoS)
    if i>=1:
      guess=guess+step
      error[i]=diff_thrt(guess,A1,A2,V1,h1,s1,f,EoS)
      if (error[i]*error[i-1])<0:
        guess=guess-0.5*step
        P2=fsolve(diff_thrt,guess,((A1,A2,V1,h1,s1,f,EoS)))
        err=diff_thrt(P2,A1,A2,V1,h1,s1,f,EoS)
        if abs(err)>1.8e-08:
          print ("Warning Throat pressure not found")
        break
    
  return P2

def diff_isenexp(P3,A3,mass2,e2,s2,f,EoS):
  s3=s2
  Prop3=Cool_Ps(P3,s3,f,EoS)
  V3=1/Prop3[1]
  h3=Prop3[3]
  u3=(2*(e2-h3))**0.5
  mass3=A3*u3/V3
  y=mass3-mass2
  return y

def diff_Shock(Z,k1,k2,k3,f,EoS):
  P4=Z[0]
  T4=Z[1]
  Prop4=Cool_PT(P4,T4,f,EoS)
  V4=1/Prop4[1]
  h4=Prop4[3]
  u4=k1*V4
  P4f=k2-u4**2.0/V4
  h4f=k3-0.5*u4**2.0
  y=np.zeros(2)
  y[0]=P4-P4f
  y[1]=h4-h4f
  return y

def fnd_stagthr(guessP1,A1,A2,ho1,so1,f,EoS):
  P1=fsolve(diff_stagthr,guessP1,(A1,A2,ho1,so1,f,EoS))
  error=diff_stagthr(P1,A1,A2,ho1,so1,f,EoS)
  if abs(error)>=1.8e-08:
    print ('from stagnation to nozzle error change the guess of the nozzle inlet',error)
  return P1


def fnd_isenexp(guess,step,A3,mass2,e2,s2,f,EoS):
  error=np.zeros(1000000)
  for i in range (1000000):
    error[i]=diff_isenexp(guess,A3,mass2,e2,s2,f,EoS)
    if i>=1:
      guess=guess+step
      error[i]=diff_isenexp(guess,A3,mass2,e2,s2,f,EoS)
      if (error[i]*error[i-1])<0:
        guess=guess-0.5*step
        P3=fsolve(diff_isenexp,guess,(A3,mass2,e2,s2,f,EoS),xtol=1.8e-015)
        err=diff_isenexp(P3,A3,mass2,e2,s2,f,EoS)
        if abs(err)>1.8e-08:
            print ("Warning isentropic expansion pressure not found")
        break
  return P3 

def pss_Shock (guess,step,k1,k2,k3,f,EoS):
  for i in range (1000):
    PT4=fsolve(diff_Shock,(guess,300),(k1,k2,k3,f,EoS))
    err=diff_Shock(PT4,k1,k2,k3,f,EoS)
    if abs(err[0])<=1.8e-08:
      if abs(err[1])<=1.8e-08:
        break
    guess=guess+step
  return PT4

def diff_Shockloc(r,P5r,A2,A5,P2,s2,mass2,e2,gss3,stp3,gss4,stp4,stp5,f,EoS):
    #gss3=P2, stp3=-0.01*P2; gss4=0.7*P1, stp4=-0.01*P1; gss5=P4, stp5=0.005*P
    print(r)
    A3=r*A2
    P3=fnd_isenexp(gss3,stp3,A3,mass2,e2,s2,f,EoS)
    s3=s2
    Prop3=Cool_Ps(P3,s3,f,EoS)
    T3=Prop3[2]
    V3=1/Prop3[1]
    h3=Prop3[3]
    a3=Prop3[6]
    u3=(2*(e2-h3))**0.5
    M3=u3/a3
    print('Upstream Normal Shock Wave achieved',M3)
    
    k1=u3/V3
    k2=P3+u3**2.0/V3
    k3=h3+u3**2.0/2
    PT4=pss_Shock(gss4,stp4,k1,k2,k3,f,EoS)
    P4=PT4[0]
    T4=PT4[1]

    Prop4=Cool_PT(P4,T4,f,EoS)
    V4=1/Prop4[1]
    h4=Prop4[3]
    s4=Prop4[4]
    a4=Prop4[6]
    u4=k1*V4
    M4=u4/a4

    mass4=A3*u3/V3
    e4=h4+u4**2.0/2
    print('Downstream Normal Shock Wave achieved',M4)

    P5=fnd_isenexp(P4,stp5,A5,mass4,e4,s4,f,EoS)
    s5=s4
    
    Prop5=Cool_Ps(P5,s5,f,EoS)
    h5=Prop5[3]
    a5=Prop5[6]
    u5=(2*(e4-h5))**0.5
    M5=u5/a5
    print('Nozzle outlet achieved',M5)
    y=P5-P5r
    return y


















