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
from Two_phase import diff_2phase_exp, sol2PhaseMP, fnd_2phase_exp, sol2Phase, Propmix
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.integrate import quad
from scipy import interpolate
from scipy.integrate import quad

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


def Moore_2phase(Pc,Xc,uc,Ath,xc,f,grid):
        
    Ac=Moore_nozzle(xc)*Ath
    xf=7
    ne=grid
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
    emom=np.zeros(ne)
    emass=np.zeros(ne)
    eener=np.zeros(ne)
    Ai=np.zeros(ne)
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
            ci[i]=propi[5]
            Yi[i]=1-Xc
            Mi[i]=ui[i]/ci[i]
            massi[i]=rhoi[i]*ui[i]*Ac
            ei[i]=hi[i]+0.5*ui[i]**2.0
            Ai[i]=Ac
    
    
            
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
            ci[i]=propi[5]
            Yi[i]=1-X
            Mi[i]=ui[i]/ci[i]
            massi[i]=rhoi[i]*ui[i]*A
            ei[i]=hi[i]+0.5*ui[i]**2.0
            Ai[i]=A
            
            I=quad(Ix,Ai[i-1],Ai[i],(Pi[i],Pi[i-1],Ai[i],Ai[i-1]))
            I=I[0]
            mom1=(Pi[i-1]+rhoi[i-1]*ui[i-1]**2.0)*Ai[i-1]+I
            mom2=(Pi[i]+rhoi[i]*ui[i]**2.0)*Ai[i]
            emom[i]=mom2-mom1    
            emass[i]=massi[i]-massi[i-1]
            eener[i]=ei[i]-ei[i-1]
            print (i)
            
        
        
        x=x+stp
        
    
    
    e=ei[ne-1]
    mass=massi[ne-1]
    
    """
    ig,rT=plt.subplots() 
    rT.plot(xi, Pi)
    rT.set(xlabel='Lenght (-)', ylabel='Temperature (K)')
    rT.grid()
    """
    return xi,Pi,rhoi,Ti,hi,si,ui,Mi,massi,ei,e,mass,ne,Yi,Ai,emom,emass,eener

def Moses_2phase(Pc,Xc,uc,Ath,xc,f,grid):
        
    Ac=Moses_nozzle(xc)*Ath
    xf=16
    ne=grid
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
    emom=np.zeros(ne)
    emass=np.zeros(ne)
    eener=np.zeros(ne)
    Ai=np.zeros(ne)
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
            ci[i]=propi[5]
            Yi[i]=1-Xc
            Mi[i]=ui[i]/ci[i]
            massi[i]=rhoi[i]*ui[i]*Ac
            ei[i]=hi[i]+0.5*ui[i]**2.0
            Ai[i]=Ac
    
    
            
        if i>0:
            A=Moses_nozzle(x)*Ath
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
            ci[i]=propi[5]
            Yi[i]=1-X
            Mi[i]=ui[i]/ci[i]
            massi[i]=rhoi[i]*ui[i]*A
            ei[i]=hi[i]+0.5*ui[i]**2.0
            Ai[i]=A
            
            I=quad(Ix,Ai[i-1],Ai[i],(Pi[i],Pi[i-1],Ai[i],Ai[i-1]))
            I=I[0]
            mom1=(Pi[i-1]+rhoi[i-1]*ui[i-1]**2.0)*Ai[i-1]+I
            mom2=(Pi[i]+rhoi[i]*ui[i]**2.0)*Ai[i]
            emom[i]=mom2-mom1    
            emass[i]=massi[i]-massi[i-1]
            eener[i]=ei[i]-ei[i-1]
            print ("STP",xi[i])
            
        
        
        x=x+stp
        
    
    
    e=ei[ne-1]
    mass=massi[ne-1]
    
    """
    ig,rT=plt.subplots() 
    rT.plot(xi, Pi)
    rT.set(xlabel='Lenght (-)', ylabel='Temperature (K)')
    rT.grid()
    """
    return xi,Pi,rhoi,Ti,hi,si,ui,Mi,massi,ei,e,mass,ne,Yi,Ai,emom,emass,eener

