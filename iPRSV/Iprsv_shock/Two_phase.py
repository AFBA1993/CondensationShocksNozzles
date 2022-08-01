# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 23:20:45 2019

@author: USER
"""

import numpy as np
from scipy.optimize import fsolve
import math
from LibiPRSV import iPRSV_mas_Ps, iPRSV_mas_PT, iPRSV_mas_PT_sat,iPRSV_mas_VT
from LibiPRSVsat import fnd_psatiPRSV, fnd_Tsat
import multiprocessing
import time


def diff_2phase_exp(Pguess,P1,X1,u1,A1,A2,f):
#   print ("Fundo",Pguess,P1,X1,u1,A1,A2,f)
    T1=fnd_Tsat(P1,f)
    propsat1=iPRSV_mas_PT_sat(f,P1,T1)
    
    rho1g=1/propsat1[2]
    h1g=propsat1[3]
    s1g=propsat1[4]
    
    
    h1l=propsat1[7]
    s1l=propsat1[8]
    
    y1=1-X1
    rho1=rho1g/X1
    h1=X1*h1g+y1*h1l
    
    s1=X1*s1g+y1*s1l
    mass1=rho1*A1*u1
    
    
    T2=fnd_Tsat(Pguess,f)
    propsat2=iPRSV_mas_PT_sat(f,Pguess,T2)
    rho2g=1/propsat2[2]
    h2g=propsat2[3]
    s2g=propsat2[4]
    
    
    h2l=propsat2[7]
    s2l=propsat2[8]
    
    X2=(s1-s2l)/(s2g-s2l)
    y2=1-X2
    h2=X2*h2g+y2*h2l
    u2=math.sqrt(2*(h1-h2+0.5*u1**2.0))
    rho2=rho2g/X2
    mass2=rho2*u2*A2
    diff=mass1-mass2
#   print("return",diff)
    return diff

def fnd_2phase_exp(Pguess,stp,u1,P1,X1,A1,A2,f):
    ne=1000000
    
    err=np.zeros(ne)
    for i in range (ne):
        err[i]=diff_2phase_exp(Pguess,P1,X1,u1,A1,A2,f)
        
        if i>1:
            if err[i]*err[i-1]<0:
                Pguess=Pguess-0.5*stp
                P2=fsolve(diff_2phase_exp,Pguess,(P1,X1,u1,A1,A2,f))
                error=diff_2phase_exp(P2,P1,X1,u1,A1,A2,f)
                if abs (error)>1.7*1e-08:
                    print("Two-phase isentropic expansion failure",error)
                break
                    
        Pguess=Pguess+stp
    return P2



def Propmix(P,X,f):
    
    Tmix=fnd_Tsat(P,f)
    Propsat=iPRSV_mas_PT_sat(f,P,Tmix)
    
    rhog=1/Propsat[2]
    hg=Propsat[3]
    sg=Propsat[4]
    cg=Propsat[12]
    
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



def sol2Phase(stp,u1,P1,X1,A1,A2,f):
    
    Prop1=Propmix(P1,X1,f)
    h1=Prop1[3]
    s1=Prop1[4]
       
    P2=fnd_2phase_exp(P1,stp,u1,P1,X1,A1,A2,f)
    T2=fnd_Tsat(P2,f)
    propsat2=iPRSV_mas_PT_sat(f,P2,T2)
    rho2g=1/propsat2[2]
    h2g=propsat2[3]
    s2g=propsat2[4]
    
    
    h2l=propsat2[7]
    s2l=propsat2[8]
    
    X2=(s1-s2l)/(s2g-s2l)
    y2=1-X2
    h2=X2*h2g+y2*h2l
    u2=math.sqrt(2*(h1-h2+0.5*u1**2.0))
    return P2,X2,u2


##################################################################
##################################################################
#Fnd 2ph
def fnd_diff2phase_P(stp,Pi,ne_P,P1,X1,u1,A1,A2,f,root,i):
#  print ("Entrada",P1,X1,u1,A1,A2,f)
    
   #NP=i
   
    err=np.zeros(ne_P)
    for i in range (ne_P):
        
        err[i]=diff_2phase_exp(Pi,P1,X1,u1,A1,A2,f)
        
        if i>1:
            if (err[i]*err[i-1])<0:
                Pguess=Pi-0.5*stp
                #print (NP,Pguess,"convergi")
                P2=fsolve(diff_2phase_exp,Pguess,(P1,X1,u1,A1,A2,f))
                root[0]=P2
                error=diff_2phase_exp(P2,P1,X1,u1,A1,A2,f)
                if abs (error)>1.71e-08:
                    print("Two-phase isentropic expansion failure",error)
                break
            if Pi<620:
                break
        
                    
        Pi+=stp
    

##################################################################



def fnd_2phase_expMP(P1,u1,X1,A1,A2,f):
    
    Pguess=P1
    Pf=620
    ne=820
        
    
    stp=(Pf-Pguess)/ne
    
    l=Pf-Pguess
    
    NPS=24
    l_P=l/NPS
    ne_P=int(ne/NPS)+3

    root=[0.0]*NPS
    processes =[]
    for i in range (NPS):
        if i==0:
            Pi=Pguess
                   
        if i>0:
            Pi+=l_P
        if i==(NPS-1):
            ne_P=int(ne/NPS)-1
        
        root[i]=multiprocessing.Array("d",2)
        process=multiprocessing.Process(target=fnd_diff2phase_P,args=(stp,Pi,ne_P,P1,X1,u1,A1,A2,f,root[i],i))
        process.start()
        processes.append(process)
        
    for process in processes:
        process.join()
    
    res_zeros=np.zeros(NPS)    
    
    for i in range (NPS):
        res_zeros[i]=root[i][0]
        
    zero=max(res_zeros)
    e=diff_2phase_exp(zero,P1,X1,u1,A1,A2,f)
    if e >1.8e-08:
        print ("Two-phase expansion not found",e)
        
    P2=zero
    #print(P2)
    return P2

def sol2PhaseMP(u1,P1,X1,A1,A2,f):
    
    
    
    Prop1=Propmix(P1,X1,f)
    h1=Prop1[3]
    s1=Prop1[4]
    P2=fnd_2phase_expMP(P1,u1,X1,A1,A2,f)
    T2=fnd_Tsat(P2,f)
    propsat2=iPRSV_mas_PT_sat(f,P2,T2)
    rho2g=1/propsat2[2]
    h2g=propsat2[3]
    s2g=propsat2[4]
    
    
    h2l=propsat2[7]
    s2l=propsat2[8]
    
    X2=(s1-s2l)/(s2g-s2l)
    y2=1-X2
    h2=X2*h2g+y2*h2l
    u2=math.sqrt(2*(h1-h2+0.5*u1**2.0))
    return P2,X2,u2





    








