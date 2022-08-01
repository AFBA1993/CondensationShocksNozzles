# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 14:05:46 2019

@author: sisea
"""
import math
import numpy as np
import CoolProp.CoolProp as CP
from scipy.optimize import fsolve
import CoolProp

def coefff (fluid):
    if fluid=="CarbonDioxide":
        A=1.98E+01
        B=7.344E-02
        C=-5.60E-05
        D=1.715E-08
            
    if fluid=="Water":
         A=32.24
         B=0.001924
         C=1.055E-05
         D=-3.596E-09
    if fluid=="R134a":
        A=16.7813
        B=286.357E-03
        C=-227.336E-06
        D=113.312E-09         
    return (A,B,C,D)


def k_1(fluid):
    if fluid =="Water":
        k1=-0.06635
    if fluid == "CarbonDioxide":
        k1=0.04285
    if fluid =='R134a':
        k1=0.0048
    return k1
    
    

def iPRSV_mas_PT(fluid,P,T):
    
    w=CP.PropsSI('acentric',fluid)
    Tc=CP.PropsSI('Tcrit',fluid)
    Pc=CP.PropsSI('Pcrit',fluid)
    m=CP.PropsSI('molar_mass',fluid)
    R=8.31462
    
    Tr=T/Tc
    
    ac=0.457235*R**2.0*Tc**2.0/Pc
    k0=0.378893+1.4897153*w-0.17131848*w**2.0+0.0196554*w**3.0
    k1=k_1(fluid)
    
    A=1.1
    B=0.25
    C=0.2
    D=1.2
    E=0.01
    
    k=k0+k1*(math.sqrt((A-D*(Tr+B))**2.0+E)+A-D*(Tr+B))*math.sqrt(Tr+C)
    alpha=(1+k*(1-math.sqrt(Tr)))**2.0
    a=ac*alpha
    b=0.077796*R*Tc/Pc
    
    A_A=a*P/(R*T)**2.0
    B_B=b*P/(R*T)
    
    
    a0=1
    a1=B_B-1
    a2=A_A-2*B_B-3*B_B**2.0
    a3=B_B**3.0+B_B**2.0-A_A*B_B
      
    coeff=[a0,a1,a2,a3]
    roots_Z=np.roots(coeff)
    roots_Z_check=np.zeros((2,3))
      
    roots_Z_check[0][0]=roots_Z[0].real
    roots_Z_check[1][0]=roots_Z[0].imag
      
    roots_Z_check[0][1]=roots_Z[1].real
    roots_Z_check[1][1]=roots_Z[1].imag
      
    roots_Z_check[0][2]=roots_Z[2].real
    roots_Z_check[1][2]=roots_Z[2].imag
     
    msg="Non two-phase zone"
    for i in range (3):
         if roots_Z_check[1][i]!=0:
            for j in range(3):
               if roots_Z_check[1][j]==0:
                  Zg=roots_Z_check[0][j]
                  Vg=R*T*Zg/P
                  break
         if roots_Z_check[1][0]==0:
            if roots_Z_check[1][1]==0:
               if  roots_Z_check[1][2]==0:
                  msg="Warning Two-phase zone"
                  Zg=max(roots_Z_check[0][0],roots_Z_check[0][1],roots_Z_check[0][2])
                  Zl=min(roots_Z_check[0][0],roots_Z_check[0][1],roots_Z_check[0][2])
                  Vg=R*T*Zg/P
                  Vl=R*T*Zl/P
    
    "Reference state"
    To=273.15+25
    Po=1e+05
    
    coef=coefff(fluid)
    A1=coef[0]
    B1=coef[1]
    C1=coef[2]
    D1=coef[3]
    
    hoideal=(A1*To+0.5*B1*To**2+1/3*C1*To**3+1/4*D1*To**4)
    hideal=(A1*T+0.5*B1*T**2+1/3*C1*T**3+1/4*D1*T**4)
    del_ideal=hideal-hoideal
    
    
    Tau_1=1-math.sqrt(Tr)
    Tau_2=math.sqrt(Tc*T)
    Tx=A-D*(Tr+B)
    Ty=math.sqrt(Tx**2.0+E)
    Tz=math.sqrt(Tr+C)
    dk_dT=k1*(Tx+Ty)/Tc*(1/(2*Tz)-D*Tz/Ty)
    da_dT=2*ac*(1+Tau_1*k)*(Tau_1*dk_dT-k/(2*Tau_2))
    
    h_r=R*T*(Zg-1)+(T*da_dT-a)/(2*math.sqrt(2)*b)*math.log((Zg+(1+math.sqrt(2))*B_B)/(Zg+(1-math.sqrt(2))*B_B))
   
    h_r1=((Zg-1)-(A_A/(B_B*8**0.5))*(1+(k*Tr**0.5)/(alpha**0.5))*math.log((Zg+(1+2**0.5)*B_B)/(Zg+(1-2**0.5)*B_B)))*R*T
    h=del_ideal+h_r
    
    polS=A1*math.log(T/To)+B1*(T-To)+C1*0.5*(T**2-To**2)+D1*1/3.0*(T**3-To**3)
    sideal=polS-R*math.log(P/Po)
    
    
    s_r=R*math.log(Zg-B_B)+(da_dT)/(2*math.sqrt(2)*b)*math.log((Zg+(1+math.sqrt(2))*B_B)/(Zg+(1-math.sqrt(2))*B_B))
    s_r1=(math.log(Zg-B_B)-(A_A/(B_B*8**0.5))*(k*Tr**0.5/alpha**0.5)*math.log((Zg+(1+2**0.5)*B_B)/(Zg+(1-2**0.5)*B_B)))*R
    
    s=s_r+sideal
    
    
    dP_dT_V=R/(Vg-b)-da_dT/(Vg*(Vg+b)+b*(Vg-b))
    dP_dV_T=-R*T/(Vg-b)**2.0+(2*a*(Vg+b))/(Vg*(Vg+b)+b*(Vg-b))**2.0
    
    #print('sr',100-s_r/s_r1*100)
    #print('hr',100-h_r/h_r1*100)
    
    
    d2k_dT2=-k1*(Tx+Ty)/Tc**2.0*(D/(Ty*Tz)+1/(4*Tz**3.0)+D**2.0*Tz/Ty**3.0*(Tx-Ty))
    d2a_dT2=2*ac*((Tau_1*dk_dT-k/(2*Tc*Tr**0.5))**2.0+(1+Tau_1*k)*(k/(4*Tc**2.0*Tr**(3/2))-dk_dT/(Tc*Tr**0.5)+Tau_1*d2k_dT2))
    Int=d2a_dT2/(2*2**0.5*b)*math.log((Zg+(1+math.sqrt(2))*B_B)/(Zg+(1-math.sqrt(2))*B_B))
    
#    d2P_dT2_V=-d2a_dT2/(Vg*(Vg+b)+b*(Vg-b))
    
    Cp_ig=A1+B1*T+C1*T**2.0+D1*T**3.0
    Cv_ig=Cp_ig-R
    Cv=Cv_ig+T*Int
    Cp=Cv-T*dP_dT_V**2.0/dP_dV_T
    
    gamma=Cp/Cv
    Ss=Vg*(-gamma*dP_dV_T)**0.5
    Ss=(Ss**2/m)**0.5
    
    "Checking speed of sound"
    a_dI=da_dT
    dA_dT_P=P/(R**2.0*T**2.0)*(a_dI-(2*a)/T)
    dB_dT_P=(-b*P)/(R*T**2.0)
    dZ_dT_P=(dA_dT_P*(B_B-Zg)+dB_dT_P*(6*B_B*Zg+2*Zg-3*B_B**2.0-2*B_B+A_A-Zg**2.0))/(3*Zg**2.0+2*(B_B-1)*Zg+(A_A-2*B_B-3*B_B**2.0))
    dV_dT_P1=R/P*(T*dZ_dT_P+Zg)
    a_dII=d2a_dT2
    Cv_r=(T*a_dII)/(b*8**0.5)*math.log((Zg+B_B*(1+2**0.5))/(Zg+B_B*(1-2**0.5)))
    Cv_1=Cv_r+Cv_ig
    Cp_r1=Cv_r+T*dP_dT_V*dV_dT_P1-R
    Cp_1=Cp_r1+Cp_ig
    Ss1=Vg*(-Cp_1/Cv_1*dP_dV_T)**0.5
    Ss1=(Ss1**2/m)**0.5
    #print('speedsound',100-Ss/Ss1*100)
    Cp=Cp/m
    V=Vg/m
    h=h/m
    s=s/m
    
    #Gibbs FREE ENERGY
    g=h-T*s
    g_r=(h_r-T*s_r)/m
    #g_ig=(del_ideal-T*sideal)/m
    
  
    return (P,V,T,h,s,Zg,Ss,msg,gamma,Cv,Cp)


def iPRSV_mas_Ps(fluid,P,s):
  def solver (T,fluid,P,s):
        si=iPRSV_mas_PT(fluid,P,T)
        si=si[4]
        y=si-s
        return y
  T=fsolve(solver,380,(fluid,P,s))
  err=solver(T,fluid,P,s)
  if abs(err)>1.8e-08:
     print('Warning PS failure: error',err)
  res=iPRSV_mas_PT(fluid,P,T)
  return res


def iPRSV_mas_hs(fluid,h,s,Pgss):
   def solver (P,fluid,h,s):
        hi=iPRSV_mas_Ps(fluid,P,s)
        hi=hi[3]
        y=hi-h
        return y
   P=fsolve(solver,Pgss,(fluid,h,s))
   err=solver(P,fluid,h,s)
   if abs(err)>1.8e-08:
      print("PR_mas_hs change P guess value")
   res=iPRSV_mas_Ps(fluid,P,s)
   return res

def iPRSV_mas_PT_sat(fluid,P,T):
    
    w=CP.PropsSI('acentric',fluid)
    Tc=CP.PropsSI('Tcrit',fluid)
    Pc=CP.PropsSI('Pcrit',fluid)
    m=CP.PropsSI('molar_mass',fluid)
    R=8.31462
    
    Tr=T/Tc
    
    ac=0.457235*R**2.0*Tc**2.0/Pc
    k0=0.378893+1.4897153*w-0.17131848*w**2.0+0.0196554*w**3.0
    k1=k_1(fluid)
    
    A=1.1
    B=0.25
    C=0.2
    D=1.2
    E=0.01
    
    k=k0+k1*(math.sqrt((A-D*(Tr+B))**2.0+E)+A-D*(Tr+B))*math.sqrt(Tr+C)
    alpha=(1+k*(1-math.sqrt(Tr)))**2.0
    a=ac*alpha
    b=0.077796*R*Tc/Pc
    
    A_A=a*P/(R*T)**2.0
    B_B=b*P/(R*T)
    
    
    a0=1
    a1=B_B-1
    a2=A_A-2*B_B-3*B_B**2.0
    a3=B_B**3.0+B_B**2.0-A_A*B_B
      
    coeff=[a0,a1,a2,a3]
    roots_Z=np.roots(coeff)
    roots_Z_check=np.zeros((2,3))
      
    roots_Z_check[0][0]=roots_Z[0].real
    roots_Z_check[1][0]=roots_Z[0].imag
      
    roots_Z_check[0][1]=roots_Z[1].real
    roots_Z_check[1][1]=roots_Z[1].imag
      
    roots_Z_check[0][2]=roots_Z[2].real
    roots_Z_check[1][2]=roots_Z[2].imag
     
    msg="Non two-phase zone"
    for i in range (3):
         if roots_Z_check[1][i]!=0:
            for j in range(3):
               if roots_Z_check[1][j]==0:
                  Zg=roots_Z_check[0][j]
                  Vg=R*T*Zg/P
                  break
         if roots_Z_check[1][0]==0:
            if roots_Z_check[1][1]==0:
               if  roots_Z_check[1][2]==0:
                  msg="Warning Two-phase zone"
                  Zg=max(roots_Z_check[0][0],roots_Z_check[0][1],roots_Z_check[0][2])
                  Zl=min(roots_Z_check[0][0],roots_Z_check[0][1],roots_Z_check[0][2])
                  Vg=R*T*Zg/P
                  Vl=R*T*Zl/P
    
    "Reference state"
    To=273.15+25
    Po=1e+05
    
    coef=coefff(fluid)
    A1=coef[0]
    B1=coef[1]
    C1=coef[2]
    D1=coef[3]
    
    hoideal=(A1*To+0.5*B1*To**2+1/3*C1*To**3+1/4*D1*To**4)
    hideal=(A1*T+0.5*B1*T**2+1/3*C1*T**3+1/4*D1*T**4)
    del_ideal=hideal-hoideal
    
    
    Tau_1=1-math.sqrt(Tr)
    Tau_2=math.sqrt(Tc*T)
    Tx=A-D*(Tr+B)
    Ty=math.sqrt(Tx**2.0+E)
    Tz=math.sqrt(Tr+C)
    dk_dT=k1*(Tx+Ty)/Tc*(1/(2*Tz)-D*Tz/Ty)
    da_dT=2*ac*(1+Tau_1*k)*(Tau_1*dk_dT-k/(2*Tau_2))
    
    h_rg=R*T*(Zg-1)+(T*da_dT-a)/(2*math.sqrt(2)*b)*math.log((Zg+(1+math.sqrt(2))*B_B)/(Zg+(1-math.sqrt(2))*B_B))
   
    
    hg=del_ideal+h_rg
    
    polS=A1*math.log(T/To)+B1*(T-To)+C1*0.5*(T**2-To**2)+D1*1/3.0*(T**3-To**3)
    sideal=polS-R*math.log(P/Po)
    
    
    s_rg=R*math.log(Zg-B_B)+(da_dT)/(2*math.sqrt(2)*b)*math.log((Zg+(1+math.sqrt(2))*B_B)/(Zg+(1-math.sqrt(2))*B_B))
    sg=s_rg+sideal
    
    h_rl=R*T*(Zl-1)+(T*da_dT-a)/(2*math.sqrt(2)*b)*math.log((Zl+(1+math.sqrt(2))*B_B)/(Zl+(1-math.sqrt(2))*B_B))
    hl=del_ideal+h_rl
    
    s_rl=R*math.log(Zl-B_B)+(da_dT)/(2*math.sqrt(2)*b)*math.log((Zl+(1+math.sqrt(2))*B_B)/(Zl+(1-math.sqrt(2))*B_B))
    sl=s_rl+sideal
    
    
    ###Speed of sound liquid
    
    dk_dT=k1*(Tx+Ty)/Tc*(1/(2*Tz)-D*Tz/Ty)
    da_dT=2*ac*(1+Tau_1*k)*(Tau_1*dk_dT-k/(2*Tau_2))
    
    
    dP_dT_V=R/(Vl-b)-da_dT/(Vl*(Vl+b)+b*(Vl-b))
    dP_dV_T=-R*T/(Vl-b)**2.0+(2*a*(Vl+b))/(Vl*(Vl+b)+b*(Vl-b))**2.0
    
    #print('sr',100-s_r/s_r1*100)
    #print('hr',100-h_r/h_r1*100)
    
    
    d2k_dT2=-k1*(Tx+Ty)/Tc**2.0*(D/(Ty*Tz)+1/(4*Tz**3.0)+D**2.0*Tz/Ty**3.0*(Tx-Ty))
    d2a_dT2=2*ac*((Tau_1*dk_dT-k/(2*Tc*Tr**0.5))**2.0+(1+Tau_1*k)*(k/(4*Tc**2.0*Tr**(3/2))-dk_dT/(Tc*Tr**0.5)+Tau_1*d2k_dT2))
    Int=d2a_dT2/(2*2**0.5*b)*math.log((Zl+(1+math.sqrt(2))*B_B)/(Zl+(1-math.sqrt(2))*B_B))

    
    Cp_ig=A1+B1*T+C1*T**2.0+D1*T**3.0
    Cv_ig=Cp_ig-R
    Cv=Cv_ig+T*Int
    Cp=Cv-T*dP_dT_V**2.0/dP_dV_T
    
    gamma=Cp/Cv
    Ssl=Vl*(-gamma*dP_dV_T)**0.5
    Ssl=(Ssl**2/m)**0.5
#    print('speedsound',100-Ss/Ss1*100)
    
    
    Vl=Vl/m
    hl=hl/m
    sl=sl/m
    Vg=Vg/m
    hg=hg/m
    sg=sg/m
    return (P,T,Vg,hg,sg,Zg,Vl,hl,sl,Zl,msg,Ssl)




def iPRSV_mas_VT(fluid,V,T):
    
    w=CP.PropsSI('acentric',fluid)
    Tc=CP.PropsSI('Tcrit',fluid)
    Pc=CP.PropsSI('Pcrit',fluid)
    m=CP.PropsSI('molar_mass',fluid)
    R=8.31462
    
    Tr=T/Tc
    
    ac=0.457235*R**2.0*Tc**2.0/Pc
    k0=0.378893+1.4897153*w-0.17131848*w**2.0+0.0196554*w**3.0
    k1=k_1(fluid)
    
    A=1.1
    B=0.25
    C=0.2
    D=1.2
    E=0.01
    
    k=k0+k1*(math.sqrt((A-D*(Tr+B))**2.0+E)+A-D*(Tr+B))*math.sqrt(Tr+C)
    alpha=(1+k*(1-math.sqrt(Tr)))**2.0
    a=ac*alpha
    b=0.077796*R*Tc/Pc
    
    V=V*m
    
    P=(R*T)/(V-b)-a/(V**2.0+2*b*V-b**2.0)
    A_A=a*P/(R*T)**2.0
    B_B=b*P/(R*T)
    Vg=V
    Zg=P*Vg/(R*T)
    
    
    
    "Reference state"
    To=273.15+25
    Po=1e+05
    
    coef=coefff(fluid)
    A1=coef[0]
    B1=coef[1]
    C1=coef[2]
    D1=coef[3]
    
    hoideal=(A1*To+0.5*B1*To**2+1/3*C1*To**3+1/4*D1*To**4)
    hideal=(A1*T+0.5*B1*T**2+1/3*C1*T**3+1/4*D1*T**4)
    del_ideal=hideal-hoideal
    
    
    Tau_1=1-math.sqrt(Tr)
    Tau_2=math.sqrt(Tc*T)
    Tx=A-D*(Tr+B)
    Ty=math.sqrt(Tx**2.0+E)
    Tz=math.sqrt(Tr+C)
    dk_dT=k1*(Tx+Ty)/Tc*(1/(2*Tz)-D*Tz/Ty)
    da_dT=2*ac*(1+Tau_1*k)*(Tau_1*dk_dT-k/(2*Tau_2))
    
    h_r=R*T*(Zg-1)+(T*da_dT-a)/(2*math.sqrt(2)*b)*math.log((Zg+(1+math.sqrt(2))*B_B)/(Zg+(1-math.sqrt(2))*B_B))
   
    h_r1=((Zg-1)-(A_A/(B_B*8**0.5))*(1+(k*Tr**0.5)/(alpha**0.5))*math.log((Zg+(1+2**0.5)*B_B)/(Zg+(1-2**0.5)*B_B)))*R*T
    h=del_ideal+h_r
    
    polS=A1*math.log(T/To)+B1*(T-To)+C1*0.5*(T**2-To**2)+D1*1/3.0*(T**3-To**3)
    sideal=polS-R*math.log(P/Po)
    
    
    s_r=R*math.log(Zg-B_B)+(da_dT)/(2*math.sqrt(2)*b)*math.log((Zg+(1+math.sqrt(2))*B_B)/(Zg+(1-math.sqrt(2))*B_B))
    s_r1=(math.log(Zg-B_B)-(A_A/(B_B*8**0.5))*(k*Tr**0.5/alpha**0.5)*math.log((Zg+(1+2**0.5)*B_B)/(Zg+(1-2**0.5)*B_B)))*R
    
    s=s_r+sideal
    
    
    dP_dT_V=R/(Vg-b)-da_dT/(Vg*(Vg+b)+b*(Vg-b))
    dP_dV_T=-R*T/(Vg-b)**2.0+(2*a*(Vg+b))/(Vg*(Vg+b)+b*(Vg-b))**2.0
    
    #print('sr',100-s_r/s_r1*100)
    #print('hr',100-h_r/h_r1*100)
    
    
    d2k_dT2=-k1*(Tx+Ty)/Tc**2.0*(D/(Ty*Tz)+1/(4*Tz**3.0)+D**2.0*Tz/Ty**3.0*(Tx-Ty))
    d2a_dT2=2*ac*((Tau_1*dk_dT-k/(2*Tc*Tr**0.5))**2.0+(1+Tau_1*k)*(k/(4*Tc**2.0*Tr**(3/2))-dk_dT/(Tc*Tr**0.5)+Tau_1*d2k_dT2))
    Int=d2a_dT2/(2*2**0.5*b)*math.log((Zg+(1+math.sqrt(2))*B_B)/(Zg+(1-math.sqrt(2))*B_B))
    
#    d2P_dT2_V=-d2a_dT2/(Vg*(Vg+b)+b*(Vg-b))
    
    Cp_ig=A1+B1*T+C1*T**2.0+D1*T**3.0
    Cv_ig=Cp_ig-R
    Cv=Cv_ig+T*Int
    Cp=Cv-T*dP_dT_V**2.0/dP_dV_T
    
    gamma=Cp/Cv
    Ss=Vg*(-gamma*dP_dV_T)**0.5
    Ss=(Ss**2/m)**0.5
    
    "Checking speed of sound"
    a_dI=da_dT
    dA_dT_P=P/(R**2.0*T**2.0)*(a_dI-(2*a)/T)
    dB_dT_P=(-b*P)/(R*T**2.0)
    dZ_dT_P=(dA_dT_P*(B_B-Zg)+dB_dT_P*(6*B_B*Zg+2*Zg-3*B_B**2.0-2*B_B+A_A-Zg**2.0))/(3*Zg**2.0+2*(B_B-1)*Zg+(A_A-2*B_B-3*B_B**2.0))
    dV_dT_P1=R/P*(T*dZ_dT_P+Zg)
    a_dII=d2a_dT2
    Cv_r=(T*a_dII)/(b*8**0.5)*math.log((Zg+B_B*(1+2**0.5))/(Zg+B_B*(1-2**0.5)))
    Cv_1=Cv_r+Cv_ig
    Cp_r1=Cv_r+T*dP_dT_V*dV_dT_P1-R
    Cp_1=Cp_r1+Cp_ig
    Ss1=Vg*(-Cp_1/Cv_1*dP_dV_T)**0.5
    Ss1=(Ss1**2/m)**0.5
    #print('speedsound',100-Ss/Ss1*100)
    
    
    V=Vg/m
    h=h/m
    s=s/m
    return (P,V,T,h,s,Zg,Ss)












"""
f='R134a'
P1=3e+06
T1=359.1
prop1=iPRSV_mas_PT(f,P1,T1)
V1=prop1[1]
s1=prop1[4]
h1=prop1[3]
ss1=prop1[6]


P2=3e+06
T2=360
prop2=iPRSV_mas_PT(f,P2,T2)
s2=prop2[4]
h2=prop2[3]
deltas=s2-s1
deltah=h2-h1



HEOS =CoolProp.AbstractState("PR",f)
HEOS.specify_phase(CoolProp.iphase_gas)
HEOS.update(CoolProp.PT_INPUTS ,P1,T1)
V1c=1/HEOS.rhomass()
h1c=HEOS.hmass()
s1c=HEOS.smass()
ss1c=HEOS.speed_sound()


HEOS.update(CoolProp.PT_INPUTS ,P2,T2)
V2c=1/HEOS.rhomass()
h2c=HEOS.hmass()
s2c=HEOS.smass()

diffs=100-(s2c-s1c)/(s2-s1)*100
print (diffs)


diffh=100-(h2c-h1c)/(h2-h1)
print (diffh)

diffvol=100-V1c/V1*100
print (diffvol)
"""


