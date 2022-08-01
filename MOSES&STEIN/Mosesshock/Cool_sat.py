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
from Coollib import Cool_Tsat, Cool_PT

def p_sat(T,f,EoS):
    HEOS = CoolProp.AbstractState(EoS, f)
    HEOS.update(CoolProp.QT_INPUTS, 1, T)
    psat=HEOS.p()
    return psat

def T_sat(P,f,EoS):
    HEOS = CoolProp.AbstractState(EoS, f)
    HEOS.update(CoolProp.QT_INPUTS, 1, T)
    psat=HEOS.T()
    return psat



def Ncl_rate2(P,T,f,EoS):
    
  Prop=Cool_PT(P,T,f,EoS)
  gamma=Prop[7]
  h_g=Prop[3]
  rho_g=Prop[1]
  s_g=Prop[4]
  cv_g=Prop[8]
  g_g=h_g-T*s_g
  
  
  
  Prop_sat=Cool_Tsat(T,f,EoS)
  Psat=Prop_sat[0]
  S=P/Psat
  
  "Properties of the liquid"
  rho_l=1/Prop_sat[6]
  h_l=Prop_sat[7]
  s_l=Prop_sat[8]
  g_l=h_l-T*s_l
  
  "Properties of the saturated vapor"
  rho_gs=1/Prop_sat[2]
  h_gs=Prop_sat[3]
  s_gs=Prop_sat[4]
  g_gs=h_gs-T*s_gs
  
  
  gs=g_g-g_gs-1/rho_l*(P-Psat)
  
  
  
  
  "Surface Tension Calculation"
  Tc=CP.PropsSI('Tcrit',f)
  
  B=235.8/1000
  b=-0.625
  miu=1.256
  
  sft=B*((Tc-T)/Tc)**miu*(1+b*(Tc-T)/Tc)
  #sft=CP.PropsSI('surface_tension','P',P,'T',T,f)
  r_crit=2*sft/(rho_l*gs)
  delG=16/3*math.pi*sft**3.0*(1/(rho_l*gs))**2.0
  
  "Nucleation rate calculation"
  m=CP.PropsSI('molar_mass',f)
  R=8.31462/m
  kb=1.3807*10**(-23)
  mh2o=2.99046*10**(-26)
  mh2o=2.991460866e-26
  qc=0.3

  Jcl=qc*rho_g/rho_l*((2*sft)/(math.pi*mh2o**3.0))**0.5*math.exp(-(delG)/(kb*T))
  #phinom=(cv_g+kb/(2*mh2o))*kb*T**2.0
  #phiden=phinom+mh2o*(T*(s_g-s_l)-(kb*T/(2*mh2o))-(2*sft/(rho_l*r_crit**2.0)))**2.0
  #phi=phinom/phiden
  C=1+2*(gamma-1)/(gamma+1)*(h_g-h_l)/(R*T)*((h_g-h_l)/(R*T)-0.5)
  J=Jcl/C

  print (J)
  return J


def Ncl_rate(P,T,f,EoS):
  "quitar 0.87 exp delG Nrate"
  Prop=Cool_PT(P,T,f,EoS)
  gamma=Prop[7]
  h_g=Prop[3]
  rho_g=Prop[1]  
  
  "Saturation properties"
  Prop_sat=Cool_Tsat(T,f,EoS)
  Psat=Prop_sat[0]
  #h_g=Prop_sat[3]
  #rho_g=1/Prop_sat[2]
  rho_l=1/Prop_sat[6]
  h_l=Prop_sat[7]
  S=P/Psat
  
  
  "Surface Tension Calculation"
  Tc=CP.PropsSI('Tcrit',f)
  
  B=235.8/1000
  b=-0.625
  miu=1.256
  
  sft=B*((Tc-T)/Tc)**miu*(1+b*(Tc-T)/Tc)
  #sft=CP.PropsSI('surface_tension','P',P,'T',T,f)
  
  "Nucleation rate calculation"
  m=CP.PropsSI('molar_mass',f)
  R=8.31462/m
  kb=1.3807*10**(-23)
  mh2o=2.99046*10**(-26)
  mh2o=2.991460866e-26
  r_crit=2*sft/(rho_l*R*T*math.log(S))

  #r_crit=2*sft/(rho_l*R*T*(math.log(S)+Gr_PT-Gr_Ps)-(P-Psat))
  
  C=1+2*(gamma-1)/(gamma+1)*(h_g-h_l)/(R*T)*((h_g-h_l)/(R*T)-0.5)
  Jcl=1*rho_g**2.0/(C*rho_l)*math.sqrt((2*sft)/(math.pi*mh2o**3))*math.exp((-4*0.87*math.pi*r_crit**2.0*sft)/(3*kb*T))
  return Jcl

def Ncl_rate3(P,T,f,EoS):
    
  Prop=Cool_PT(P,T,f,EoS)
  gamma=Prop[7]
  h_g=Prop[3]
  rho_g=Prop[1]
  s_g=Prop[4]
  cv_g=Prop[8]
  g_g=h_g-T*s_g
  
  
  
  Prop_sat=Cool_Tsat(T,f,EoS)
  Psat=Prop_sat[0]
  S=P/Psat
  
  "Properties of the liquid"
  rho_l=1/Prop_sat[6]
  h_l=Prop_sat[7]
  s_l=Prop_sat[8]
  g_l=h_l-T*s_l
  
  "Properties of the saturated vapor"
  rho_gs=1/Prop_sat[2]
  h_gs=Prop_sat[3]
  s_gs=Prop_sat[4]
  g_gs=h_gs-T*s_gs
  
  
  gs=g_g-g_gs-1/rho_l*(P-Psat)
  
  
  
  
  "Surface Tension Calculation"
  Tc=CP.PropsSI('Tcrit',f)
  
  B=235.8/1000
  b=-0.625
  miu=1.256
  
  sft=B*((Tc-T)/Tc)**miu*(1+b*(Tc-T)/Tc)
  #sft=CP.PropsSI('surface_tension','P',P,'T',T,f)
  r_crit=2*sft/(rho_l*gs)
  delG=16/3*math.pi*sft**3.0*(1/(rho_l*gs))**2.0
  
  "Nucleation rate calculation"
  m=CP.PropsSI('molar_mass',f)
  R=8.31462/m
  kb=1.3807*10**(-23)
  mh2o=2.99046*10**(-26)
  mh2o=2.991460866e-26
  qc=1

  Jcl=qc*rho_g/rho_l*((2*sft)/(math.pi*mh2o**3.0))**0.5*math.exp(-(delG)/(kb*T))
  phinom=(cv_g+kb/(2*mh2o))*kb*T**2.0
  phiden=phinom+mh2o*(T*(s_g-s_l)-(kb*T/(2*mh2o)))**2.0
  #-(2*sft/(rho_l*r_crit**2.0
  phi=phinom/phiden
  #C=1+2*(gamma-1)/(gamma+1)*(h_g-h_l)/(R*T)*((h_g-h_l)/(R*T)-0.5)
  #J=Jcl/C
  J=rho_gs/rho_g*phi*Jcl

  print ("eu")
  return J ,r_crit