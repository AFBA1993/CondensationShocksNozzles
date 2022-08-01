# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 18:28:54 2019

@author: USER
"""
import CoolProp
from scipy.optimize import fsolve, newton
import numpy as np

def Cool_PT(P,T,f,EoS):
    #print(T)
    #if T<200:
       #T=218
    
    HEOS = CoolProp.AbstractState(EoS, f)
    HEOS.specify_phase(CoolProp.iphase_gas)
    HEOS.update(CoolProp.PT_INPUTS, P, T)
    
    rho=HEOS.rhomass()
    h=HEOS.hmass()
    s=HEOS.smass()
    cp=HEOS.cpmass()
    cv=HEOS.cvmass()
    R=HEOS.gas_constant()/HEOS.molar_mass()
    Z=P/(rho*R*T)
    gamma=cp/cv
    c=HEOS.speed_sound()
    return P,rho,T,h,s,Z,c,gamma,cv
"""
def Cool_Ps(P,s,f,EoS,gssT,lim):
   def solver (T,P,s,f,EoS,gssT,lim):
        #lim=494
        if T<lim:
            T=lim
            
        #print ("yo ",P,T)
        si=Cool_PT(P,T,f,EoS)
        si=si[4]
        y=si-s
        return y
    
   guessT=0.985*gssT
   T=fsolve(solver,guessT,(P,s,f,EoS,gssT,lim))
   err=solver(T,P,s,f,EoS,gssT,lim)
   if abs(err)>1.8e-08:
     print('Warning PS failure: error',err, P)
   res=Cool_PT(P,T,f,EoS)
   return res
"""

def Cool_hs(h,s,Pgss,f,EoS,uplimT):
   def solver (P,h,s,F,EoS,uplimT):
        hi=Cool_Ps(P,s,f,EoS,uplimT)
        hi=hi[3]
        y=hi-h
        return y
   P=fsolve(solver,Pgss,(h,s,f,EoS,uplimT))
   err=solver(P,h,s,f,EoS,uplimT)
   if abs(err)>1.8e-08:
      print("PR_mas_hs change P guess value",err)
   res=Cool_Ps(P,s,f,EoS,uplimT)
   return res

def Cool_Tsat(T,f,EoS):
    HEOS = CoolProp.AbstractState(EoS, f)
    HEOS.update(CoolProp.QT_INPUTS, 1, T)
    P=HEOS.p()
    Vg=1/HEOS.rhomass()
    hg=HEOS.hmass()
    sg=HEOS.smass()
    R=HEOS.gas_constant()/HEOS.molar_mass()
    Zg=P/(1/Vg*R*T)
    
    HEOS.update(CoolProp.QT_INPUTS, 0, T)
    Vl=1/HEOS.rhomass()
    hl=HEOS.hmass()
    sl=HEOS.smass()
    Zl=P/(1/Vg*R*T)
 
    return (P,T,Vg,hg,sg,Zg,Vl,hl,sl,Zl)
    
def Cool_Psat(P,f,EoS):
    HEOS = CoolProp.AbstractState(EoS, f)
    HEOS.update(CoolProp.PQ_INPUTS, P, 1)
    T=HEOS.T()
    Vg=1/HEOS.rhomass()
    hg=HEOS.hmass()
    sg=HEOS.smass()
    cg=HEOS.speed_sound()
    R=HEOS.gas_constant()/HEOS.molar_mass()
    Zg=P/(1/Vg*R*T)
    
    HEOS.update(CoolProp.PQ_INPUTS, P, 0)
    Vl=1/HEOS.rhomass()
    hl=HEOS.hmass()
    sl=HEOS.smass()
    cl=HEOS.speed_sound()
    Zl=P/(1/Vg*R*T)
    
    return (P,T,Vg,hg,sg,Zg,Vl,hl,sl,Zl,cg,cl)



def Cool_Ps_imp(P,s,f,EoS,upT,dT):
   def solver (T,P,s,f,EoS):
        print ("yo ",T)
        si=Cool_PT(P,T,f,EoS)
        si=si[4]
        y=si-s
        return y
    
def Cool_Ps(P,s,f,EoS,uplimT):
    
     def solver (T,P,s,f,EoS):
        si=Cool_PT(P,T,f,EoS)
        si=si[4]
        y=si-s
        return y
     T=uplimT
     stp=0.5
     err=np.zeros(1000)
     for i in range(1000):
         err[i]=solver(T,P,s,f,EoS)
         if i>0:
             if err[i]*err[i-1]<0:
                 gssT=T-0.5*stp
                 res=fsolve(solver,gssT,(P,s,f,EoS))
                 e=solver(res,P,s,f,EoS)
                 if abs(e)>1.8e-07:
                     print ("Ps failure",err)
                 break
         
         
         T-=stp
     res=Cool_PT(P,res,f,EoS)
     return res    
    





"""
f='Water'
EoS='HEOS'


To=273.15+366
Po=101*100000

propo=Cool_PT(Po,To,f,EoS)
so=propo[4]

s1=so
P1=0.5*Po



pro1o=Cool_Ps(P1,s1,f,EoS,To)


s1=pro1o[4]
err=so-s1
print (err)
"""



















