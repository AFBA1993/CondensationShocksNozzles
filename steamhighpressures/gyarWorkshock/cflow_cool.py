# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 22:06:51 2019

@author: USER
"""
import numpy as np
from scipy.optimize import fsolve
from Coollib import Cool_PT, Cool_Ps



def diff_thrt(P2,A1,A2,V1,h1,s1,f,EoS,gssT,lim):
    s2=s1
    prop2=Cool_Ps(P2,s2,f,EoS,gssT,lim)
    V2=1/prop2[1]
    h2=prop2[3]
    a2=prop2[6]
    e2=h2+a2**2.0/2
    mass2=A2*a2/V2
    u1=(2*(e2-h1))**0.5
    mass1=A1*u1/V1
    y=mass1-mass2
    return y

def diff_stagthr(P1,A1,A2,ho1,so1,f,EoS,gssT,lim):
  s1=so1
  Propo1=Cool_Ps(P1,s1,f,EoS,gssT,lim)
  V1=1/Propo1[1]
  h1=Propo1[3]
  s1=Propo1[4]
  
  step=-0.001*P1
  guess=0.8*P1
  P2=fnd_throat(guess,step,A1,A2,V1,h1,s1,f,EoS,gssT,lim)
  s2=s1
  
  Prop2=Cool_Ps(P2,s2,f,EoS,gssT,lim)
  
  a2=Prop2[6]
  h2=Prop2[3]
  e2=h2+a2**2.0/2
  u1=(2*(e2-h1))**0.5
  e1=h1+u1**2.0/2
  y=ho1-e1
  return y


def fnd_throat(guess,step,A1,A2,V1,h1,s1,f,EoS,gssT,lim):
  error=np.zeros(1000)
  for i in range (1000):
    #print(guess)
    error[i]=diff_thrt(guess,A1,A2,V1,h1,s1,f,EoS,gssT,lim)
    if i>=1:
      guess=guess+step
      error[i]=diff_thrt(guess,A1,A2,V1,h1,s1,f,EoS,gssT,lim)
      if (error[i]*error[i-1])<0:
        guess=guess-0.5*step
        P2=fsolve(diff_thrt,guess,((A1,A2,V1,h1,s1,f,EoS,gssT,lim)))
        err=diff_thrt(P2,A1,A2,V1,h1,s1,f,EoS,gssT,lim)
        if abs(err)>1.8e-08:
          print ("Warning Throat pressure not found",err)
        break
    
  return P2

def diff_isenexp(P3,A3,mass2,e2,s2,f,EoS,gssT,lim):
  s3=s2
  Prop3=Cool_Ps(P3,s3,f,EoS,gssT,lim)
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

def fnd_stagthr(guessP1,A1,A2,ho1,so1,f,EoS,gssT,lim):
  P1=fsolve(diff_stagthr,guessP1,(A1,A2,ho1,so1,f,EoS,gssT,lim))
  error=diff_stagthr(P1,A1,A2,ho1,so1,f,EoS,gssT,lim)
  if abs(error)>=1.8e-08:
    print ('from stagnation to nozzle error change the guess of the nozzle inlet',error)
  return P1


def fnd_isenexp(guess,step,A3,mass2,e2,s2,f,EoS,gssT,lim):
  error=np.zeros(1000000)
  for i in range (1000000):
    error[i]=diff_isenexp(guess,A3,mass2,e2,s2,f,EoS,gssT,lim)
    if i>=1:
      guess=guess+step
      error[i]=diff_isenexp(guess,A3,mass2,e2,s2,f,EoS,gssT,lim)
      if (error[i]*error[i-1])<0:
        guess=guess-0.5*step
        P3=fsolve(diff_isenexp,guess,(A3,mass2,e2,s2,f,EoS,gssT,lim),xtol=1e-015)
        err=diff_isenexp(P3,A3,mass2,e2,s2,f,EoS,gssT,lim)
        if abs(err)>1.8e-08:
            print ("Warning isentropic expansion pressure not found",err)
        break
  return P3 




"""
f="Water"
EoS="HEOS"
A1=B4M(-20)*1
A2=1

To=273.15+366
Po=101*100000
propo=Cool_PT(Po,To,f,EoS)
so=propo[4]
ho=propo[3]

P1=0.98*Po
s1=so
prop1=Cool_Ps(P1,so,f,EoS)
T1=prop1[2]
V1=1/prop1[1]
h1=prop1[3]

P2=fnd_throat(0.8*P1,-0.001*P1,A1,A2,V1,h1,s1,f,EoS)
s2=s1
  
Prop2=Cool_Ps(P2,s2,f,EoS)
  
a2=Prop2[6]
h2=Prop2[3]
e2=h2+a2**2.0/2
u1=(2*(e2-h1))**0.5
e1=h1+u1**2.0/2
y=ho-e1
print(y)


P1x=fnd_stagthr(Po,A1,A2,ho,so,f,EoS)
"""

