# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 19:59:56 2019

@author: USER
"""

import numpy as np
from Coollib import Cool_Psat
from scipy.optimize import fsolve
import math


def diff_2phase_exp(Pguess,P1,X1,u1,A1,A2,f,EoS):
    
    Prop1=Cool_Psat(P1,f,EoS)
    rho1g=1/Prop1[2]
    h1g=Prop1[3]
    s1g=Prop1[4]
    
    h1l=Prop1[7]
    s1l=Prop1[8]
    
    y1=1-X1
    rho1=rho1g/X1
    h1=X1*h1g+y1*h1l
    
    s1=X1*s1g+y1*s1l
    mass1=rho1*A1*u1
    
    
    Prop2=Cool_Psat(Pguess,f,EoS)
    rho2g=1/Prop2[2]
    h2g=Prop2[3]
    s2g=Prop2[4]
    
    h2l=Prop2[7]
    s2l=Prop2[8]
    
    X2=(s1-s2l)/(s2g-s2l)
    y2=1-X2
    h2=X2*h2g+y2*h2l
    u2=math.sqrt(2*(h1-h2+0.5*u1**2.0))
    rho2=rho2g/X2
    mass2=rho2*u2*A2
    
    diff=mass1-mass2
    return diff

def fnd_2phase_exp(Pguess,stp,u1,P1,X1,A1,A2,f,EoS):
    ne=1000000
    
    err=np.zeros(ne)
    for i in range (ne):
        err[i]=diff_2phase_exp(Pguess,P1,X1,u1,A1,A2,f,EoS)
        if i>1:
            if err[i]*err[i-1]<0:
                Pguess=Pguess-0.5*stp
                P2=fsolve(diff_2phase_exp,Pguess,(P1,X1,u1,A1,A2,f,EoS))
                error=diff_2phase_exp(P2,P1,X1,u1,A1,A2,f,EoS)
                if abs (error)>1.7*1e-08:
                    print("Two-phase isentropic expansion failure",error)
                break
                    
        Pguess=Pguess+stp
    return P2

def Propmix(P,X,f,EoS):
    Propsat=Cool_Psat(P,f,EoS)
    Tmix=Propsat[1]
    
    rhog=1/Propsat[2]
    hg=Propsat[3]
    sg=Propsat[4]
    cg=Propsat[10]
    
    #rhol=1/Propsat[6]
    hl=Propsat[7]
    sl=Propsat[8]
    cl=Propsat[11]
    
    Pmix=P
    rhomix=rhog/X
    hmix=hg*X+(1-X)*hl
    smix=sg*X+(1-X)*sl
    cmix=cg*X+(1-X)*cl
    return Pmix, rhomix, Tmix, hmix, smix, cmix

def sol2Phase(stp,u1,P1,X1,A1,A2,f,EoS):
    
    Prop1=Propmix(P1,X1,f,EoS)
    h1=Prop1[3]
    s1=Prop1[4]
       
    P2=fnd_2phase_exp(P1,stp,u1,P1,X1,A1,A2,f,EoS)
    Prop2=Cool_Psat(P2,f,EoS)

    h2g=Prop2[3]
    s2g=Prop2[4]
    
    h2l=Prop2[7]
    s2l=Prop2[8]

    
    X2=(s1-s2l)/(s2g-s2l)
    y2=1-X2
    h2=X2*h2g+y2*h2l
    u2=math.sqrt(2*(h1-h2+0.5*u1**2.0))
    return P2,X2,u2














