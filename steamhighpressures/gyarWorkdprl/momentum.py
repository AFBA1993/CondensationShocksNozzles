#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 20:57:51 2020

@author: andres777
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import curve_fit

def integrand(x,coeff):
    
    y=coeff[0]*x**20+coeff[1]*x**19+coeff[2]*x**18+coeff[3]*x**17+coeff[4]*x**16+coeff[5]*x**15+coeff[6]*x**14+coeff[7]*x**13+coeff[8]*x**12+coeff[9]*x**11+coeff[10]*x**10+coeff[11]*x**9+coeff[12]*x**8+coeff[13]*x**7+coeff[14]*x**6+coeff[15]*x**5+coeff[16]*x**4+coeff[17]*x**3+coeff[18]*x**2+coeff[19]*x+coeff[20]
    #y=coeff[0]*x**11+coeff[1]*x**10+coeff[2]*x**9+coeff[3]*x**8+coeff[4]*x**7+coeff[5]*x**6+coeff[6]*x**5+coeff[7]*x**4+coeff[8]*x**3+coeff[9]*x**2+coeff[10]*x+coeff[11]
    return y

def err_mom(ne,stp,ii,results):
    "SETTING LIMITS"
    Pi=results["P"][ii]
    ui=results["u"][ii]
    rhoi=results["rho"][ii]
    Ai=results["A"][ii]
    
    
    
    Pf=results["P"][ii+stp]
    uf=results["u"][ii+stp]
    rhof=results["rho"][ii+stp]
    Af=results["A"][ii+stp]
    
   
    AA=np.zeros(ne)
    PP=np.zeros(ne)
    
    kk=ii
    for j in range(ne):
        AA[j]=results["A"][kk]
        PP[j]=results["P"][kk]
        kk+=1
        
    
    pf=np.polyfit(AA,PP,ne)
    xx=np.linspace(Ai,Af,ne)
    yy=np.polyval(pf,xx)
    
    gp=20
    
    coeffn=np.zeros(gp+1)
    j=ne
    for l in range (ne+1):
        coeffn[gp]=pf[j]
        j=j-1
        gp=gp-1
    
    
    I=quad(integrand,Ai,Af,args=(coeffn))
    I=I[0]+I[1]
    
    err_mom=Pf*Af+rhof*uf**2.0*Af-(Pi*Ai+rhoi*ui**2.0*Ai+I)
    return err_mom

def counter(l,xcond,results):
    for i in range(100):
        xs=results["x"][l]
        if xs==xcond:
            break
        l+=1
    return i
            
   

def Err_mom(xcond,xff,results):
    
    
    ne=4
    stp=ne-1
    
    ii=0
    
    
    #errm=np.zeros(10000)
    errTm=np.zeros(10000)
    xf_i=np.zeros(10000)
    c=0
    
    for i in range (10000):
        xi=results["x"][ii]
        xf=results["x"][ii+stp]
        
        if xf>=xcond:
            ne_f=counter(ii,xcond,results)
            stp_f=ne_f-1
            
            #errm[i]=err_mom(ne_f,stp_f,ii)
            errTm[c]=err_mom(ne_f,stp_f,ii,results)
            xf=results["x"][ii+stp_f]
            xf_i[c]=xf
            break        
            
        #errm[i]=err_mom(ne,stp,ii)
        errTm[c]=err_mom(ne,stp,ii,results)
        xf_i[c]=xf
        ii+=stp
        c+=1
    ii=ii+ne_f+1
    ne=2
    stp=ne-1
    
    #errm2=np.zeros(10000)
    for i in range (10000):
        xi=results["x"][ii]
        xf=results["x"][ii+stp]
        
        if xf>=xff:
            #errm2[i]=err_mom(ne,stp,ii)
            errTm[c]=err_mom(ne,stp,ii,results)
            xf_i[c]=xf
            break        
            
        #errm2[i]=err_mom(ne,stp,ii)
        errTm[c]=err_mom(ne,stp,ii,results)
        xf_i[c]=xf
        ii+=stp
        c+=1
    
    
    sol_Emom=np.zeros(c)
    sol_xf=np.zeros(c)
    
    for i in range (c):
        sol_Emom[i]=errTm[i]
        sol_xf[i]=xf_i[i]
    
    errmax=max(sol_Emom)
    errmin=min(sol_Emom)




    plt.plot(sol_xf,sol_Emom)
    return sol_Emom,sol_xf,errmax,errmin








