# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 22:06:51 2019

@author: USER
"""
import numpy as np
from scipy.optimize import fsolve
from Coollib import Cool_PT, Cool_Ps
import multiprocessing 
import time 


def diff_thrt(P2,A1,A2,V1,h1,s1,f,EoS,uplimT):
    s2=s1
    prop2=Cool_Ps(P2,s2,f,EoS,uplimT)
    V2=1/prop2[1]
    h2=prop2[3]
    a2=prop2[6]
    e2=h2+a2**2.0/2
    #print(e2-h1)
    mass2=A2*a2/V2
    #print (mass2)
    u1=(2*(e2-h1))**0.5
    #print (u1)
    mass1=A1*u1/V1
    y=mass1-mass2
    j=y.imag

    #print(y,P2)
    
    return y

def diff_stagthr(P1,A1,A2,ho1,so1,f,EoS,uplimT):
  s1=so1
  Propo1=Cool_Ps(P1,s1,f,EoS,uplimT)
  V1=1/Propo1[1]
  h1=Propo1[3]
  s1=Propo1[4]
  
  step=-0.00005*P1
  guess=0.7*P1
  #P2=fnd_throat(guess,step,A1,A2,V1,h1,s1,f,EoS,uplimT)
  P2=fnd_throat_MP(P1,A1,A2,V1,h1,s1,f,EoS,uplimT)
  s2=s1
  
  Prop2=Cool_Ps(P2,s2,f,EoS,uplimT)
  
  a2=Prop2[6]
  h2=Prop2[3]
  e2=h2+a2**2.0/2
  u1=(2*(e2-h1))**0.5
  e1=h1+u1**2.0/2
  y=ho1-e1
  return y


def fnd_throat(guess,step,A1,A2,V1,h1,s1,f,EoS,uplimT):
  error=np.zeros(100000)
  for i in range (100000):
    #print(guess)
    error[i]=diff_thrt(guess,A1,A2,V1,h1,s1,f,EoS,uplimT)
    if i>0:
      if (error[i]*error[i-1])<0:
        guess=guess-0.5*step
        P2=fsolve(diff_thrt,guess,((A1,A2,V1,h1,s1,f,EoS,uplimT)))
        err=diff_thrt(P2,A1,A2,V1,h1,s1,f,EoS,uplimT)
        print("garganta enconrada")
        
        if abs(err)>1.8e-08:
          print ("Warning Throat pressure not found",err)
        break
    #print (guess)
    guess=guess+step
    
  return P2

def diff_isenexp(P3,A3,mass2,e2,s2,f,EoS,uplimT):
  s3=s2
  Prop3=Cool_Ps(P3,s3,f,EoS,uplimT)
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

def fnd_stagthr(guessP1,A1,A2,ho1,so1,f,EoS,uplimT):
    
  
  P1=fsolve(diff_stagthr,guessP1,(A1,A2,ho1,so1,f,EoS,uplimT))
  error=diff_stagthr(P1,A1,A2,ho1,so1,f,EoS,uplimT)
  if abs(error)>=1.8e-08:
    print ('from stagnation to nozzle error change the guess of the nozzle inlet',error)
  return P1


def fnd_isenexp(guess,step,A3,mass2,e2,s2,f,EoS,uplimT):
  error=np.zeros(1000000)
  for i in range (1000000):
    error[i]=diff_isenexp(guess,A3,mass2,e2,s2,f,EoS,uplimT)
    if i>=1:
      guess=guess+step
      error[i]=diff_isenexp(guess,A3,mass2,e2,s2,f,EoS,uplimT)
      if (error[i]*error[i-1])<0:
        guess=guess-0.5*step
        P3=fsolve(diff_isenexp,guess,(A3,mass2,e2,s2,f,EoS,uplimT))
        err=diff_isenexp(P3,A3,mass2,e2,s2,f,EoS,uplimT)
        if abs(err)>1.8e-08:
            print ("Warning isentropic expansion pressure not found",err)
        break
  return P3 


########################################################################################
################################################ MUTIPROCESSING
#####Process function 
  
def fnd_thrt_pf(Pi,stp,ne_P,A1,A2,V1,h1,s1,f,EoS,uplimT,root,i):
    
    err=np.zeros(ne_P)
    for i in range (ne_P):
        #Pi=float(Pi)
        err[i]=diff_thrt(Pi,A1,A2,V1,h1,s1,f,EoS,uplimT)
        if i>0:
            if (err[i]*err[i-1])<0:
                gssP=Pi-0.5*stp
                root[0]=fsolve(diff_thrt,gssP,(A1,A2,V1,h1,s1,f,EoS,uplimT))
                
                break
        Pi+=stp
        #print(c,Pi,err[i])



#######################################################################################



def fnd_throat_MP(P1,A1,A2,V1,h1,s1,f,EoS,uplimT): 
    Pi=0.7*P1
    Pf=0.52*P1
    LT=Pf-Pi
    ne=2800
    stp=LT/ne
    #val=3506580.0+1
    #j=diff_thrt(val,A1,A2,V1,h1,s1,f,EoS,To)
    
    
    NPS=4
    
    ne_P=int(ne/NPS)+2
    
    process=[]
    root=[0.0]*NPS
    
    for i in range (NPS):
        if i==0:
            Pi=Pi
        if i>0:
            Pi+=LT/NPS
        if i==(NPS-1):
            ne_P=int(ne/NPS)
        
        root[i]=multiprocessing.Array("d",2)
        p=multiprocessing.Process(target=fnd_thrt_pf,args=(Pi,stp,ne_P,A1,A2,V1,h1,s1,f,EoS,uplimT,root[i],i))
        p.start()
        process.append(p)
        
    for p in process:
        p.join()    
    
        
    res_zeros=[0.0]*NPS
    for i in range (NPS):
        res_zeros[i]=root[i][0]   
        
    P2=max(res_zeros)
    err=diff_thrt(P2,A1,A2,V1,h1,s1,f,EoS,uplimT)
    print("Thrt found for",P1,err)

    
    return (P2)
 
"""
f="CarbonDioxide"
EoS="HEOS"

P1=0.99*6440000.000000001
s1=1810.9905644993078
A1=4.115058537377536
A2=1
To=313

t1=time.time()

Prop1=Cool_Ps(P1,s1,f,EoS,To)
V1=1/Prop1[1]
h1=Prop1[3]
P2=fnd_throat_MP(P1,A1,A2,V1,h1,s1,f,EoS,To)
"""






