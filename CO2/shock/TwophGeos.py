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
from scipy import interpolate

def Ix(x,P4,P3,A4,A3):
        m=(P4-P3)/(A4-A3)
        b=P3-m*A3
        y=m*x+b       
        return y
"""
def B2_Bier(x):
    if x<=2.01355789988531:
        y=0.008123539812761996*x**3.0+0.012445354905944056*x**2.0+0.009423166697206373*x+1.0000000000000002
    if x>2.01355789988531:
        y=-0.0010506010852403273*x**3.0+0.02646759857633423*x**2.0+0.025729978042863963*x+0.9852091434987051

    return y 

"""
def B2_Bier(x):
    x_data=np.zeros(4)
    y_data=np.zeros(4)
    
    x_data[0]=0    
    x_data[1]=2.25352643886249
    x_data[2]=4.28656558799125
    x_data[3]=6
    
    y_data[0]=1
    y_data[1]=1.16005380805109
    y_data[2]=1.4973850292927
    y_data[3]=1.86585703150537
    
    tck = interpolate.splrep(x_data, y_data)
    return interpolate.splev(x, tck)

def B2_Bier_2phase(Pc,Xc,uc,Ath,xc,f,EoS):
        
    Ac=B2_Bier(xc)*Ath
    xf=6
    ne=400
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
            A=B2_Bier(x)*Ath
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
            A=B2_Bier(x)*Ath
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
            
        
        print(xi[i])
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














































def Moore_2phase(Pc,Xc,uc,Ath,xc,f,EoS):
        
    Ac=Moore_nozzle(xc)*Ath
    xf=7
    ne=100
    diff=xf-xc
    stp=diff/ne
    stpg=-0.005*Pc #always is the problem
    

    
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
            
        if i>0:
            A=Moore_nozzle(x)*Ath
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
            print (i)
        
        
        x=x+stp
        
    
    
    e=ei[ne-1]
    mass=massi[ne-1]
    
    ig,rT=plt.subplots() 
    rT.plot(xi, Mi)
    rT.set(xlabel='Lenght (-)', ylabel='Temperature (K)')
    rT.grid() 
    return xi,Pi,rhoi,Ti,hi,si,ui,Mi,massi,ei,e,mass,ne,Yi 











"""
f="Water"
EoS="HEOS"
Pc=2.836150737845585900e+04
Xc=9.720612697089644483e-01
uc=5.678796273617973611e+02
Ath=1
xc=75.0

exp=Barschdorff_2phase(Pc,Xc,uc,Ath,xc,f,EoS)

"""









    
