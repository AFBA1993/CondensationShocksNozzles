# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 20:04:14 2019

@author: USER
"""
import CoolProp
import CoolProp.CoolProp as CP
from scipy.optimize import fsolve
import numpy as np
import math
from Coollib import Cool_Tsat, Cool_PT, Cool_Psat

def p_sat(T,f,EoS):
    HEOS = CoolProp.AbstractState(EoS, f)
    HEOS.update(CoolProp.QT_INPUTS, 1, T)
    psat=HEOS.p()
    return psat

def T_sat(P,f,EoS):
    HEOS = CoolProp.AbstractState(EoS, f)
    HEOS.update(CoolProp.PQ_INPUTS, P, 1)
    Tsat=HEOS.T()
    return Tsat





def Ncl_rate(P,T,f,EoS,k):
  Prop=Cool_PT(P,T,f,EoS)
  gamma=Prop[7]
  h_g=Prop[3]
  rho_g=Prop[1]  
  
  "Saturation properties"
  Prop_sat=Cool_Tsat(T,f,EoS)
  Psat=Prop_sat[0]
  
  Prop_sat=Cool_Psat(P,f,EoS)
  rho_l=1/Prop_sat[6]
  h_l=Prop_sat[7]
  S=P/Psat
  
  
  "Surface Tension Calculation"
  Tc=CP.PropsSI('Tcrit',f)
  
  B=235.8/1000
  b=-0.625
  miu=1.256
  
  sft=B*((Tc-T)/Tc)**miu*(1+b*(Tc-T)/Tc)
    
  "Nucleation rate calculation"
  m=CP.PropsSI('molar_mass',f)
  R=8.31462/m
  kb=1.3807*10**(-23)
  mh2o=2.99046*10**(-26)
  mh2o=2.991460866e-26
  r_crit=2*sft/(rho_l*R*T*math.log(S))

    
  C=1+2*(gamma-1)/(gamma+1)*(h_g-h_l)/(R*T)*((h_g-h_l)/(R*T)-0.5)
  Jcl=1*rho_g**2.0/(C*rho_l)*math.sqrt((2*sft)/(math.pi*mh2o**3))*math.exp((k*-4*math.pi*r_crit**2.0*sft)/(3*kb*T))
  
  
  return Jcl, r_crit




