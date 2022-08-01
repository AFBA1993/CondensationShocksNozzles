# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 19:13:33 2019

@author: USER
"""
import numpy as np
from scipy.optimize import fsolve
from Coollib import Cool_PT, Cool_Ps, Cool_hs, Cool_Tsat, Cool_Psat
from cflow_cool import  fnd_isenexp, fnd_throat,pss_Shock, diff_Shockloc,diff_isenexp
import matplotlib.pyplot as plt
from Cool_sat import p_sat, Ncl_rate, T_sat
from scipy.integrate import quad
from scipy import interpolate

def f1(x):
    x_points = [ 2, 2.3399999999999998579, 2.8020000000000000462, 3.3000000000000002665, 3.7400000000000002132, 4.3419999999999996376, 5.3390000000000004121]
    y_points = [4,3.404,2.75,2.188,1.784,1.406,1.086]

    tck = interpolate.splrep(x_points, y_points,k=2)
    return interpolate.splev(x, tck)
def f2(x):
    x_points = [ 5.3390000000000004121,8.22]
    y_points = [1.086, 1]

    tck = interpolate.splrep(x_points, y_points,k=1)
    return interpolate.splev(x, tck)
def f3(x):
    x_points = [8.22, 9.34, 10.188, 11.104, 11.96, 13.398,16]
    y_points = [1,1.024, 1.064, 1.128, 1.22,1.414,1.94186]

    tck = interpolate.splrep(x_points, y_points)
    return interpolate.splev(x, tck)

def Moses_nozzle(x):
    if (x>=2) & (x<=5.3390000000000004121):
        y=f1(x)
    if (x>5.3390000000000004121) & (x<=8.22):
        y=f2(x)
    if (x>8.22):
        y=f3(x)
    return y
"""
def Moses_nozzle(x):
    #x0=3.3000000000000003 has 8.22
    coeff=np.zeros(9)
    coeff[0]=7.693190168029795e-06
    coeff[1]=-0.00027245053954241366
    coeff[2]=0.003410535151559491
    coeff[3]=-0.01095737299571049
    coeff[4]=-0.14476226891860317
    coeff[5]=1.6754496274573842
    coeff[6]=-7.042439005537108
    coeff[7]=12.204146601389098
    coeff[8]=-3.4999366366261913
    
    y=coeff[0]*x**8.0+coeff[1]*x**7.0+coeff[2]*x**6.0+coeff[3]*x**5.0+coeff[4]*x**4.0+coeff[5]*x**3.0+coeff[6]*x**2.0+coeff[7]*x+coeff[8]
    
    if x>8.22:
        y=0.000631275589*x**3.0+-0.005897299805*x**2.0-0.02120488737*x+1.222157092


    

    return y
"""    



def Ix(x,P4,P3,A4,A3):
        m=(P4-P3)/(A4-A3)
        b=P3-m*A3
        y=m*x+b       
        return y



def Moses_expansion(step,A2,P2,s2,f,EoS,Jcl,xg,k):
  prop2=Cool_Ps(P2,s2,f,EoS)
  V2=1/prop2[1]
  h2=prop2[3]
  a2=prop2[6]
  u2=a2
  e2=h2+a2**2.0/2
  mass2=A2*a2/V2
  xth=8.22
  L=16-6
  div=62*10
  stp=1/div
  ne=int(div*L)+2
     
  Cool_rho=np.zeros(ne)
  Cool_x=np.zeros(ne)
  Cool_r=np.zeros(ne)
  Cool_P=np.zeros(ne)
  Cool_V=np.zeros(ne)
  Cool_T=np.zeros(ne)
  Cool_h=np.zeros(ne)
  Cool_s=np.zeros(ne)
  Cool_u=np.zeros(ne)
  Cool_M=np.zeros(ne)
  Cool_A=np.zeros(ne)
  Cool_S=np.zeros(ne)
  Cool_Jcl=np.zeros(ne)
  Cool_Sub=np.zeros(ne)
  Cool_rcrit=np.zeros(ne)
  ei=np.zeros(ne)
  massi=np.zeros(ne)
  emass=np.zeros(ne)
  eener=np.zeros(ne)
  emom=np.zeros(ne)
  x=6
   
  for i in range (ne):
    r=Moses_nozzle(x)
    Ai=r*A2
    if i==0:
      Pi=fnd_isenexp(P2,0.07*P2,Ai,mass2,e2,s2,f,EoS)
      si=s2
      propi=Cool_Ps(Pi,si,f,EoS)
      ui=(2*(h2+0.5*u2**2-propi[3]))**0.5 
    if x==xth:
      r=1
      Pi=P2
      si=s2
      propi=Cool_Ps(Pi,si,f,EoS)
      ui=propi[6]
    if i>0:
      Pi=fnd_isenexp(Cool_P[i-1],-0.0001*Cool_P[i-1],Ai,mass2,e2,s2,f,EoS)
      si=s2
      propi=Cool_Ps(Pi,si,f,EoS)
      ui=(2*(h2+0.5*u2**2-propi[3]))**0.5
    
      
    if x>L:
      break
  
      
    Mi=ui/propi[6]
    Cool_A[i]=Ai
    Cool_x[i]=x
    Cool_r[i]=r
    Cool_P[i]=propi[0]
    Cool_V[i]=1/propi[1]
    Cool_rho[i]=propi[1]
    Cool_T[i]=propi[2]
    Cool_h[i]=propi[3]
    Cool_s[i]=propi[4]
    Cool_u[i]=ui
    Cool_M[i]=Mi
    ei[i]=Cool_h[i]+0.5*Cool_u[i]**2.0
    massi[i]=Cool_rho[i]*Cool_u[i]*Cool_A[i]
    
    if i>0:
        I=I=quad(Ix,Cool_A[i-1],Cool_A[i],(Cool_P[i],Cool_P[i-1],Cool_A[i],Cool_A[i-1]))
        I=I[0]
        mom1=(Cool_P[i-1]+Cool_rho[i-1]*Cool_u[i-1]**2.0)*Cool_A[i-1]+I
        mom2=(Cool_P[i]+Cool_rho[i]*Cool_u[i]**2.0)*Cool_A[i]
        emom[i]=mom2-mom1    
        emass[i]=massi[i]-massi[i-1]
        eener[i]=ei[i]-ei[i-1]
    
    Cool_Sub[i]=T_sat(Cool_P[i],f,EoS)-Cool_T[i]
    
    
    Cool_S[i]=Cool_P[i]/p_sat(Cool_T[i],f,EoS)
    if Cool_S[i]>1:
        Ji=Ncl_rate(Cool_P[i],Cool_T[i],f,EoS,k)
        Cool_Jcl[i]=Ji[0]
        rcrit=Ji[1]
        Cool_rcrit=Ji[1]
    
    
    print(Cool_x[i],Cool_S[i],Cool_Jcl[i])
    
    x_cond=Cool_x[i]
    u_cond=Cool_u[i]
    T_cond=Cool_T[i]
    P_cond=Cool_P[i]
    J_cond=Cool_Jcl[i]
    rcrit_cond=0
    
    if f=="Water":
      if Cool_T[i]<=273.16:       
         break
     
    if Cool_x[i]>=xg:
        break
    if Cool_Jcl[i]>0:
        break
         
   
    x=x+stp

  return Cool_x,Cool_r,Cool_P,Cool_rho,Cool_T,Cool_h,Cool_s,Cool_u,Cool_M,Cool_S,Cool_Jcl,x_cond,u_cond,T_cond,P_cond,ne,Cool_A,J_cond,rcrit_cond,Cool_Sub,Cool_rcrit,emom, emass, eener,ei




