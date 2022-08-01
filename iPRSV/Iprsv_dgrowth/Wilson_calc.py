# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 17:49:12 2019

@author: sisea
"""
import numpy as np
from scipy.optimize import fsolve
import math
from LibiPRSV import iPRSV_mas_Ps, iPRSV_mas_PT_sat,iPRSV_mas_PT
from cflowiPRSV import fnd_isenexp
import matplotlib.pyplot as plt
from LibiPRSVsat import fnd_psatiPRSV
import CoolProp.CoolProp as CP

def Ncl_rate(f,P,T,Psat):
  Prop=iPRSV_mas_PT(f,P,T)
  gamma=Prop[8]
  S=P/Psat
  
  "Saturation properties"
  Prop_sat=iPRSV_mas_PT_sat(f,Psat,T)
  rho_l=1/Prop_sat[6]
  rho_g=1/Prop_sat[2]
  
  h_l=Prop_sat[7]
  h_g=Prop_sat[3]
  
  
  "Surface Tension Calculation"
  Tc=647.15
  
  B=235.8/1000
  b=-0.625
  miu=1.256
  
  sft=B*((Tc-T)/Tc)**miu*(1+b*(Tc-T)/Tc)
  
  "Nucleation rate calculation"
  m=CP.PropsSI('molar_mass',f)
  R=8.31462/m
  kb=1.3807*10**(-23)
  mh2o=2.99046*10**(-26)
  
  r_crit=2*sft/(rho_l*R*T*math.log(S))
  C=1+2*(gamma-1)/(gamma+1)*(h_g-h_l)/(R*T)*((h_g-h_l)/(R*T)-0.5)
  Jcl=rho_g**2.0/(C*rho_l)*math.sqrt((2*sft)/(math.pi*mh2o**3))*math.exp((-4*math.pi*r_crit**2.0*sft)/(3*kb*T))
  return Jcl
















