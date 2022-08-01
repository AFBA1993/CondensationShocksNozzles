# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 21:26:52 2019

@author: USER
"""
import numpy as np
from scipy.optimize import fsolve
import math
from LibiPRSV import iPRSV_mas_Ps, iPRSV_mas_PT, iPRSV_mas_PT_sat,iPRSV_mas_VT
from LibiPRSVsat import fnd_psatiPRSV, fnd_Tsat
from Two_phase import diff_2phase_exp, fnd_2phase_exp, sol2Phase, Propmix
import matplotlib.pyplot as plt

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


def Moore_2phase(Pc,Xc,uc,Ath,xc,f):
        
    Ac=Moore_nozzle(xc)*Ath
    xf=7
    ne=100
    diff=xf-xc
    stp=diff/ne
    stpg=-0.008*Pc #always is the problem
    

    
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
    
    x=xc
    
    for i in range (ne):
        if i==0:
            xi[i]=xc
            Pi[i]=Pc
            Xi[i]=Xc
            ui[i]=uc
            propi=Propmix(Pc,Xc,f)
            rhoi[i]=propi[1]
            Ti[i]=propi[2]
            hi[i]=propi[3]
            si[i]=propi[4]
            #ci[i]=propi[5]
            Yi[i]=1-Xc
            #Mi[i]=ui[i]/ci[i]
            massi[i]=rhoi[i]*ui[i]*Ac
            ei[i]=hi[i]+0.5*ui[i]**2.0
            
        if i>0:
            A=Moore_nozzle(x)*Ath
            sol=sol2Phase(stpg,uc,Pc,Xc,Ac,A,f)
            P=sol[0]
            X=sol[1]
            u=sol[2]
            propi=Propmix(P,X,f)
            xi[i]=x
            Pi[i]=P
            Xi[i]=X
            ui[i]=u
            rhoi[i]=propi[1]
            Ti[i]=propi[2]
            hi[i]=propi[3]
            si[i]=propi[4]
            #ci[i]=propi[5]
            Yi[i]=1-X
            #Mi[i]=ui[i]/ci[i]
            massi[i]=rhoi[i]*ui[i]*A
            ei[i]=hi[i]+0.5*ui[i]**2.0
            print (i)
        
        
        x=x+stp
        
    
    
    e=ei[ne-1]
    mass=massi[ne-1]
    
    ig,rT=plt.subplots() 
    rT.plot(xi, Pi)
    rT.set(xlabel='Lenght (-)', ylabel='Temperature (K)')
    rT.grid() 
    return xi,Pi,rhoi,Ti,hi,si,ui,Mi,massi,ei,e,mass,ne,Yi



def Barschdorff_2phase(Pc,Xc,uc,Ath,xc,f):
        
    Ac=Barschdorff_nozzle(xc)*Ath
    xf=120
    ne=100
    diff=xf-xc
    stp=diff/ne
    stpg=-0.008*Pc#always is the problem
    
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
    
    x=xc
    
    for i in range (ne):
        if i==0:
            xi[i]=xc
            Pi[i]=Pc
            Xi[i]=Xc
            ui[i]=uc
            propi=Propmix(Pc,Xc,f)
            rhoi[i]=propi[1]
            Ti[i]=propi[2]
            hi[i]=propi[3]
            si[i]=propi[4]
            #ci[i]=propi[5]
            Yi[i]=1-Xc
            #Mi[i]=ui[i]/ci[i]
            massi[i]=rhoi[i]*ui[i]*Ac
            ei[i]=hi[i]+0.5*ui[i]**2.0
            
        if i>0:
            A=Barschdorff_nozzle(x)*Ath
            
            sol=sol2Phase(stpg,uc,Pc,Xc,Ac,A,f)
            P=sol[0]
            X=sol[1]
            u=sol[2]
            propi=Propmix(P,X,f)
            xi[i]=x
            Pi[i]=P
            Xi[i]=X
            ui[i]=u
            rhoi[i]=propi[1]
            Ti[i]=propi[2]
            hi[i]=propi[3]
            si[i]=propi[4]
            #ci[i]=propi[5]
            Yi[i]=1-X
            #Mi[i]=ui[i]/ci[i]
            massi[i]=rhoi[i]*ui[i]*A
            ei[i]=hi[i]+0.5*ui[i]**2.0
        
        print (x,A)
        x=x+stp
        
    e=ei[ne-1]
    mass=massi[ne-1]    
    
    ig,rT=plt.subplots() 
    rT.plot(xi, Ti)
    rT.set(xlabel='Lenght (-)', ylabel='Temperature (K)')
    rT.grid() 
    return xi,Pi,rhoi,Ti,hi,si,ui,Mi,massi,ei,e,mass,ne,Yi



"""
f="Water"
Pc=25000
Xc=0.993
uc=485
Ath=1
xc=56

l=Barschdorff_2phase(Pc,Xc,uc,Ath,xc,f)
"""