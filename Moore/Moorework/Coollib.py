# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 18:28:54 2019

@author: USER
"""
import CoolProp
from scipy.optimize import fsolve


def Cool_DT(rho,T,f,EoS):
    HEOS = CoolProp.AbstractState(EoS, f)
    HEOS.specify_phase(CoolProp.iphase_gas)
    HEOS.update(CoolProp.DmassT_INPUTS, rho, T)
    P=HEOS.p()
    rho=HEOS.rhomass()
    h=HEOS.hmass()
    s=HEOS.smass()
    cp=HEOS.cpmass()
    cv=HEOS.cvmass()
    R=HEOS.gas_constant()/HEOS.molar_mass()
    
    gamma=cp/cv
    c=HEOS.speed_sound()
    lamda=HEOS.conductivity()
    mu=HEOS.viscosity()
    R=HEOS.gas_constant()/HEOS.molar_mass()
    Z=P/(rho*R*T)
    return P,rho,T,h,s,Z,c,gamma,cv,cp,lamda,mu,R
    
    

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
    lamda=HEOS.conductivity()
    mu=HEOS.viscosity()
    R=HEOS.gas_constant()/HEOS.molar_mass()
    return P,rho,T,h,s,Z,c,gamma,cv,cp,lamda,mu,R

def Cool_Ps(P,s,f,EoS):
   def solver (T,P,s,f,EoS):
        si=Cool_PT(P,T,f,EoS)
        si=si[4]
        y=si-s
        return y
   T=fsolve(solver,320,(P,s,f,EoS))
   err=solver(T,P,s,f,EoS)
   if abs(err)>1.8e-08:
     print('Warning PS failure: error',err)
   res=Cool_PT(P,T,f,EoS)
   return res

def Cool_hs(h,s,Pgss,f,EoS):
   def solver (P,h,s,F,EoS):
        hi=Cool_Ps(P,s,f,EoS)
        hi=hi[3]
        y=hi-h
        return y
   P=fsolve(solver,Pgss,(h,s,f,EoS))
   err=solver(P,h,s,f,EoS)
   if abs(err)>1.8e-08:
      print("PR_mas_hs change P guess value",err)
   res=Cool_Ps(P,s,f,EoS)
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
    cg=HEOS.speed_sound()
    
    HEOS.update(CoolProp.QT_INPUTS, 0, T)
    Vl=1/HEOS.rhomass()
    hl=HEOS.hmass()
    sl=HEOS.smass()
    Zl=P/(1/Vg*R*T)
    cl=HEOS.speed_sound()
 
    return (P,T,Vg,hg,sg,Zg,Vl,hl,sl,Zl,cg,cl)
    
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
    cpl=HEOS.cpmass()
    
    
    return (P,T,Vg,hg,sg,Zg,Vl,hl,sl,Zl,cg,cl,cpl)





