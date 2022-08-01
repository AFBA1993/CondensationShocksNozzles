# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 00:31:28 2019

@author: USER
"""
import numpy as np
from Coollib import Cool_Psat, Cool_PT
from scipy.optimize import fsolve
from Twophase import fnd_2phase_exp, Propmix, sol2Phase
import math
import matplotlib.pyplot as plt
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



def Moore_nozzle(x):
  xth=2
  if x<=xth:
    Ar=-0.0625*x+1.127
  if x>xth:
    Ar=0.088*x+0.824
  return Ar

def Barschdorff_nozzle(x):
  A=3.00E-15
  B=1.00E-12
  C=2.00E-10
  D=8.00E-09
  E=3.00E-05
  F=2.00E-06
  G=1.00E+00
  xth=0
  if x<xth:
    y=A*x**6.0+B*x**5.0+C*x**4.0+D*x**3.0+E*x**2+F*x+G
  
  if x>=xth:
    x=-x
    y=A*x**6.0+B*x**5.0+C*x**4.0+D*x**3.0+E*x**2+F*x+G
  
  return y



def Barschdorff_2phase(Pc,Xc,uc,Ath,xc,f,EoS):
        
    Ac=Barschdorff_nozzle(xc)*Ath
    xf=120
    ne=100
    diff=xf-xc
    stp=diff/ne
    stpg=-0.0001*Pc#always is the problem
    
    xi=np.zeros(ne)
    Pi=np.zeros(ne)
    rhoi=np.zeros(ne)
    Ti=np.zeros(ne)
    hi=np.zeros(ne)
    si=np.zeros(ne)
    Xi=np.zeros(ne)
    Yi=np.zeros(ne)
    ui=np.zeros(ne)
    ci=np.zeros(ne)
    massi=np.zeros(ne)
    ei=np.zeros(ne)
    Mi=np.zeros(ne)
    Ai=np.zeros(ne)
    
    
    x=xc
    
    for i in range (ne):
        if i==0:
            A=Barschdorff_nozzle(xc)*Ath
            xi[i]=xc
            Pi[i]=Pc
            Xi[i]=Xc
            ui[i]=uc
            propi=Propmix(Pc,Xc,f,EoS)
            rhoi[i]=propi[1]
            Ti[i]=propi[2]
            hi[i]=propi[3]
            si[i]=propi[4]
            ci[i]=propi[5]
            Yi[i]=1-Xc
            Mi[i]=ui[i]/ci[i]
            massi[i]=rhoi[i]*ui[i]*Ac
            ei[i]=hi[i]+0.5*ui[i]**2.0
            Ai[i]=A
        if i>0:
            A=Barschdorff_nozzle(x)*Ath
            sol=sol2Phase(stpg,uc,Pc,Xc,Ac,A,f,EoS)
            P=sol[0]
            X=sol[1]
            u=sol[2]
            propi=Propmix(P,X,f,EoS)
            xi[i]=x
            Pi[i]=P
            Xi[i]=X
            ui[i]=u
            rhoi[i]=propi[1]
            Ti[i]=propi[2]
            hi[i]=propi[3]
            si[i]=propi[4]
            ci[i]=propi[5]
            Yi[i]=1-X
            Mi[i]=ui[i]/ci[i]
            massi[i]=rhoi[i]*ui[i]*A
            ei[i]=hi[i]+0.5*ui[i]**2.0
            Ai[i]=A
            print(xi[i])
        
        
        x=x+stp
        
    e=ei[ne-1]
    mass=massi[ne-1]    
    
    ig,rT=plt.subplots() 
    rT.plot(xi, Ti)
    rT.set(xlabel='Lenght (-)', ylabel='Temperature (K)')
    rT.grid() 
    return xi,Pi,rhoi,Ti,hi,si,ui,Mi,massi,ei,e,mass,ne,Yi,Ai

def B2M_2phase(Pc,Xc,uc,Ath,xc,f,EoS):
        
    Ac=B2M(xc)*Ath
    xf=100
    ne=50
    diff=xf-xc
    stp=diff/ne
    stpg=-0.0001*Pc#always is the problem
    
    xi=np.zeros(ne)
    Pi=np.zeros(ne)
    rhoi=np.zeros(ne)
    Ti=np.zeros(ne)
    hi=np.zeros(ne)
    si=np.zeros(ne)
    Xi=np.zeros(ne)
    Yi=np.zeros(ne)
    ui=np.zeros(ne)
    ci=np.zeros(ne)
    massi=np.zeros(ne)
    ei=np.zeros(ne)
    Mi=np.zeros(ne)
    Ai=np.zeros(ne)
    
    
    x=xc
    
    for i in range (ne):
        if i==0:
            A=B2M(xc)*Ath
            xi[i]=xc
            Pi[i]=Pc
            Xi[i]=Xc
            ui[i]=uc
            propi=Propmix(Pc,Xc,f,EoS)
            rhoi[i]=propi[1]
            Ti[i]=propi[2]
            hi[i]=propi[3]
            si[i]=propi[4]
            ci[i]=propi[5]
            Yi[i]=1-Xc
            Mi[i]=ui[i]/ci[i]
            massi[i]=rhoi[i]*ui[i]*Ac
            ei[i]=hi[i]+0.5*ui[i]**2.0
            Ai[i]=A
        if i>0:
            A=B2M(x)*Ath
            sol=sol2Phase(stpg,uc,Pc,Xc,Ac,A,f,EoS)
            P=sol[0]
            X=sol[1]
            u=sol[2]
            propi=Propmix(P,X,f,EoS)
            xi[i]=x
            Pi[i]=P
            Xi[i]=X
            ui[i]=u
            rhoi[i]=propi[1]
            Ti[i]=propi[2]
            hi[i]=propi[3]
            si[i]=propi[4]
            ci[i]=propi[5]
            Yi[i]=1-X
            Mi[i]=ui[i]/ci[i]
            massi[i]=rhoi[i]*ui[i]*A
            ei[i]=hi[i]+0.5*ui[i]**2.0
            Ai[i]=A
            print(xi[i],si[i],"area",Ai[i])
        
        
        x=x+stp
        
    e=ei[ne-1]
    mass=massi[ne-1]    
    
    ig,rT=plt.subplots() 
    rT.plot(xi, Ti)
    rT.set(xlabel='Lenght (-)', ylabel='Temperature (K)')
    rT.grid() 
    return xi,Pi,rhoi,Ti,hi,si,ui,Mi,massi,ei,e,mass,ne,Yi,Ai

def B4M_2phase(Pc,Xc,uc,Ath,xc,f,EoS):
        
    Ac=B4M(xc)*Ath
    xf=30
    ne=500
    diff=xf-xc
    stp=diff/ne
    stpg=-0.001*Pc#always is the problem
    
    xi=np.zeros(ne)
    Pi=np.zeros(ne)
    rhoi=np.zeros(ne)
    Ti=np.zeros(ne)
    hi=np.zeros(ne)
    si=np.zeros(ne)
    Xi=np.zeros(ne)
    Yi=np.zeros(ne)
    ui=np.zeros(ne)
    ci=np.zeros(ne)
    massi=np.zeros(ne)
    ei=np.zeros(ne)
    Mi=np.zeros(ne)
    Ai=np.zeros(ne)
    
    emom=np.zeros(ne)
    emass=np.zeros(ne)
    eener=np.zeros(ne)
    
    
    x=xc
    
    for i in range (ne):
        if i==0:
            A=B4M(xc)*Ath
            xi[i]=xc
            Pi[i]=Pc
            Xi[i]=Xc
            ui[i]=uc
            propi=Propmix(Pc,Xc,f,EoS)
            rhoi[i]=propi[1]
            Ti[i]=propi[2]
            hi[i]=propi[3]
            si[i]=propi[4]
            ci[i]=propi[5]
            Yi[i]=1-Xc
            Mi[i]=ui[i]/ci[i]
            massi[i]=rhoi[i]*ui[i]*Ac
            ei[i]=hi[i]+0.5*ui[i]**2.0
            Ai[i]=A
        if i>0:
            A=B4M(x)*Ath
            sol=sol2Phase(stpg,uc,Pc,Xc,Ac,A,f,EoS)
            P=sol[0]
            X=sol[1]
            u=sol[2]
            propi=Propmix(P,X,f,EoS)
            xi[i]=x
            Pi[i]=P
            Xi[i]=X
            ui[i]=u
            rhoi[i]=propi[1]
            Ti[i]=propi[2]
            hi[i]=propi[3]
            si[i]=propi[4]
            ci[i]=propi[5]
            Yi[i]=1-X
            Mi[i]=ui[i]/ci[i]
            massi[i]=rhoi[i]*ui[i]*A
            ei[i]=hi[i]+0.5*ui[i]**2.0
            Ai[i]=A
            
            I=I=quad(Ix,Ai[i-1],Ai[i],(Pi[i],Pi[i-1],Ai[i],Ai[i-1]))
            I=I[0]
            mom1=(Pi[i-1]+rhoi[i-1]*ui[i-1]**2.0)*Ai[i-1]+I
            mom2=(Pi[i]+rhoi[i]*ui[i]**2.0)*Ai[i]
            emom[i]=mom2-mom1    
            emass[i]=massi[i]-massi[i-1]
            eener[i]=ei[i]-ei[i-1]
            
        
        print(xi[i],"area",emom[i])
        x=x+stp
        
    e=ei[ne-1]
    mass=massi[ne-1]    
    """
    ig,rT=plt.subplots() 
    rT.plot(xi, Ti)
    rT.set(xlabel='Lenght (-)', ylabel='Temperature (K)')
    rT.grid() 
    """
    return xi,Pi,rhoi,Ti,hi,si,ui,Mi,massi,ei,e,mass,ne,Yi,Ai,emom,emass,eener

"""
xi=0
xf=30

l=xf-xi
ne=200
stp=l/ne
xx=np.zeros(ne)
yy=np.zeros(ne)
x=0
for i in range(ne):
    yy[i]=B4M(x)
    xx[i]=x
    
    x+=stp

plt.plot(xx,yy)
"""


































































    
