# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 21:20:33 2019

@author: USER
"""
import math
import numpy as np
import CoolProp.CoolProp as CP
from scipy.optimize import fsolve
from LibiPRSV import iPRSV_mas_PT,iPRSV_mas_PT_sat 
import multiprocessing, threading
import time


def k_1(fluid):
    if fluid =="Water":
        k1=-0.06635
    if fluid == "CarbonDioxide":
        k1=0.04285
    if fluid =='R134a':
        k1=0.0048
    return k1



def Reynoldscolonna(P,T,f):       
    w=CP.PropsSI('acentric',f)
    Tc=CP.PropsSI('Tcrit',f)
    Pc=CP.PropsSI('Pcrit',f)
    R=8.31462
    
    Tr=T/Tc
    
    ac=0.457235*R**2.0*Tc**2.0/Pc
    k0=0.378893+1.4897153*w-0.17131848*w**2.0+0.0196554*w**3.0
    k1=k_1(f)
    
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
                  break
         if roots_Z_check[1][0]==0:
            if roots_Z_check[1][1]==0:
               if  roots_Z_check[1][2]==0:
                  msg="Warning Two-phase zone"
                  Zg=max(roots_Z_check[0][0],roots_Z_check[0][1],roots_Z_check[0][2])
                  Zl=min(roots_Z_check[0][0],roots_Z_check[0][1],roots_Z_check[0][2])

                  
                  
    ln_fpl=(Zl-1)-math.log(Zl-B_B)-a/(2*2**0.5*b*R*T)*math.log((Zl+(1+2**0.5)*B_B)/(Zl+(1-2**0.5)*B_B))
    ln_fpg=(Zg-1)-math.log(Zg-B_B)-a/(2*2**0.5*b*R*T)*math.log((Zg+(1+2**0.5)*B_B)/(Zg+(1-2**0.5)*B_B))
    
    f_pl=math.exp(ln_fpl)
    f_pg=math.exp(ln_fpg)
    y=f_pg-f_pl
    return y
  
def Reynoldscolonna2(T,P,f):       
    w=CP.PropsSI('acentric',f)
    Tc=CP.PropsSI('Tcrit',f)
    Pc=CP.PropsSI('Pcrit',f)
    R=8.31462
    
    Tr=T/Tc
    
    ac=0.457235*R**2.0*Tc**2.0/Pc
    k0=0.378893+1.4897153*w-0.17131848*w**2.0+0.0196554*w**3.0
    k1=k_1(f)
    
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
                  break
         if roots_Z_check[1][0]==0:
            if roots_Z_check[1][1]==0:
               if  roots_Z_check[1][2]==0:
                  msg="Warning Two-phase zone"
                  Zg=max(roots_Z_check[0][0],roots_Z_check[0][1],roots_Z_check[0][2])
                  Zl=min(roots_Z_check[0][0],roots_Z_check[0][1],roots_Z_check[0][2])

                  
                  
    ln_fpl=(Zl-1)-math.log(Zl-B_B)-a/(2*2**0.5*b*R*T)*math.log((Zl+(1+2**0.5)*B_B)/(Zl+(1-2**0.5)*B_B))
    ln_fpg=(Zg-1)-math.log(Zg-B_B)-a/(2*2**0.5*b*R*T)*math.log((Zg+(1+2**0.5)*B_B)/(Zg+(1-2**0.5)*B_B))
    
    f_pl=math.exp(ln_fpl)
    f_pg=math.exp(ln_fpg)
    y=f_pg-f_pl
    return y


def PsatiPRSV(f,T,P,stp):
    
  
  Tc=CP.PropsSI('Tcrit',f)
  Tr=T/Tc
  if abs(Tr)>=0.99:
  
    Psat=fsolve(Reynoldscolonna,P,(T,f))
    err=Reynoldscolonna (Psat,T,f)
    if abs(err)>1.8e-08:
        print ("Saturation Pressure not found",err)
  if abs(Tr)<0.99:
      
      if abs (Tr)>0.56:
        guess=P
        error=np.zeros(100000)
        for i in range(100000):      
          error[i]=Reynoldscolonna(guess,T,f)
          if i>=1:
            guess=guess+stp
            error[i]=Reynoldscolonna (guess,T,f)
            if (error[i]*error[i-1])<0:
              guess=guess-0.5*stp
              Psat=fsolve(Reynoldscolonna,guess,(T,f))
              err=Reynoldscolonna (Psat,T,f)
              if abs(err)>1.8e-08:
                  print ("Saturation Pressure not found",err)
              break
  if abs(Tr)<=0.56:
    guess=P
    error=np.zeros(100000)
    stp=stp/96
    for i in range(100000):      
      error[i]=Reynoldscolonna(guess,T,f)
      if i>=1:
        guess=guess+stp
        error[i]=Reynoldscolonna (guess,T,f)
        if (error[i]*error[i-1])<0:
          guess=guess-0.5*stp
          Psat=fsolve(Reynoldscolonna,guess,(T,f))
          err=Reynoldscolonna (Psat,T,f)
          if abs(err)>1.8e-08:
              print ("Saturation Pressure not found",err)
          break

      

  return Psat




def lines_iPRSVsat(f):
  
  Pc=CP.PropsSI('Pcrit',f)
  Tc=CP.PropsSI('Tcrit',f)
  propC=iPRSV_mas_PT(f,Pc,Tc)
   
  P=Pc
  T=Tc
  stp=-0.5e+05
    
  ne=20
  psat=np.zeros(ne)
  Tsat=np.zeros(ne)
  Vsat_l=np.zeros(ne)
  Vsat_lr=np.zeros(ne)
  Vsat_g=np.zeros(ne)
  Vsat_gr=np.zeros(ne)
  psat_r=np.zeros(ne)
  Tsat_r=np.zeros(ne)
  s_g=np.zeros(ne)
  s_gr=np.zeros(ne)
  s_l=np.zeros(ne)
  s_lr=np.zeros(ne)
    
  for j in range (ne):
  
    if j==0:
      propC=iPRSV_mas_PT(f,P,T)
      psat[j]=propC[0]
      psat_r[j]=psat[j]/Pc
      Vc=propC[1]
      sc=propC[4]
      s_g[j]=propC[4]
      s_gr[j]=s_g[j]/sc
      
      s_l[j]=propC[4]
      s_lr[j]=s_l[j]/sc
      
      Vsat_g[j]=Vc
      Vsat_gr[j]=Vc/Vc
      
      Vsat_l[j]=Vc
      Vsat_lr[j]=Vc/Vc
  
      Tsat[j]=T
      Tsat_r[j]=Tsat[j]/Tc
      
    if j>=1:    
      for i in range (10000000):
        msg=iPRSV_mas_PT(f,P,T)
        msg=msg[7]

        #print(Pr,Tr)
        if msg=="Warning Two-phase zone":
          Psat=PsatiPRSV(f,T,P,stp)
          break
        if P<=0:
          break
        P=P+stp
      
      psat[j]=Psat
      Tsat[j]=T
      Tsat_r[j]=Tsat[j]/Tc
      propsat=iPRSV_mas_PT_sat(f,Psat,T)
      print(propsat[0]/Pc)
      
      Vsat_l[j]=propsat[6]
      Vsat_lr[j]=Vsat_l[j]/Vc
      s_g[j]=propsat[4]
      s_gr[j]=s_g[j]/sc
      
      s_l[j]=propsat[8]
      s_lr[j]=s_l[j]/sc
      
      
      Vsat_g[j]=propsat[2]
      Vsat_gr[j]=Vsat_g[j]/Vc
      
      psat_r[j]=psat[j]/Pc
      Tsat_r[j]=Tsat[j]/Tc
   
   
    T=T-1
    print(T)
  x1_vr=np.zeros(2*ne)
  y1_Pr=np.zeros(2*ne)
  
  x1_v=np.zeros(2*ne)
  y1_P=np.zeros(2*ne)
  
  x1_s=np.zeros(2*ne)
  y1_T=np.zeros(2*ne)
  
  "Vr vs Psat" 
  "liquido"   
  Vsat_reslr=sorted(Vsat_lr)
  Psat_reslr=sorted(psat_r)
  
  Vsat_resgr=Vsat_gr
  Psat_resgr=psat_r
      
  for i in range(2*ne):
    if i<=(ne-1):
      x1_vr[i]=Vsat_reslr[i]
      y1_Pr[i]=Psat_reslr[i]
    if i>(ne-1):
      x1_vr[i]=Vsat_resgr[i-ne]
      y1_Pr[i]=Psat_resgr[i-ne]
  
    
  "V vs Psat"
  Vsat_resl=sorted(Vsat_l)
  Psat_resl=sorted(psat)
  
  Vsat_resg=Vsat_g
  Psat_resg=psat
  
  for i in range(2*ne):
    if i<=(ne-1):
      x1_v[i]=Vsat_resl[i]
      y1_P[i]=Psat_resl[i]
    if i>(ne-1):
      x1_v[i]=Vsat_resg[i-ne]
      y1_P[i]=Psat_resg[i-ne]
      
  "s vs T"
  ssat_resl=sorted(s_l)
  Tsat_resl=sorted(Tsat)
  
  ssat_resg=s_g
  Tsat_resg=Tsat        
          
  for i in range(2*ne):
    if i<=(ne-1):
      x1_s[i]=ssat_resl[i]
      y1_T[i]=Tsat_resl[i]
    if i>(ne-1):
      x1_s[i]=ssat_resg[i-ne]
      y1_T[i]=Tsat_resg[i-ne]
      
  
  
  return x1_vr,y1_Pr,x1_v,y1_P,x1_s,y1_T 

def Ncl_rate(f,P,T,Psat,Tsat,k):
  Prop=iPRSV_mas_PT(f,P,T)
  gamma=Prop[8]
  rho_g=1/Prop[1]
  h_g=Prop[3]
  
  
  S=P/Psat
  
  "Saturation properties"
  Prop_sat=iPRSV_mas_PT_sat(f,P,Tsat)
  rho_l=1/Prop_sat[6]  
  h_l=Prop_sat[7]

  
  
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
  mh2o=2.991460866e-26

  r_crit=2*sft/(rho_l*R*T*math.log(S))  
  C=1+2*(gamma-1)/(gamma+1)*(h_g-h_l)/(R*T)*((h_g-h_l)/(R*T)-0.5)
  
  "Volumetric"
  Jcl=rho_g**2.0/(C*rho_l)*math.sqrt((2*sft)/(math.pi*mh2o**3))*math.exp(k*(-4*math.pi*r_crit**2.0*sft)/(3*kb*T))
  return Jcl, r_crit



def fnd_Tsat(P,f):  
  guessT=273.16
  ne=400
  err=np.zeros(ne)
  
  for i in range (ne):
    err[i]=Reynoldscolonna2(guessT,P,f)
    if i>0:
      if err[i]*err[i-1]<0:
        guessT=guessT-1
        Tsat=fsolve(Reynoldscolonna2,guessT,(P,f))
        if (Reynoldscolonna2(Tsat,P,f))>1.8e-08:
          print ('Tsat not found', (Reynoldscolonna2(Tsat,P,f)))
        break
    guessT=guessT+2
    
  return Tsat



##################################################################
def fnd_psatiPRSV(f,T):
  #trocar guess de pressao 
    P=0.1*CP.PropsSI('Pcrit',f)
    #stp=-0.6e+05
    stp=-200000
    for i in range(10000000):
        msg=iPRSV_mas_PT(f,P,T)
        msg=msg[7]
        if msg=='Warning Two-phase zone':
            Psat=PsatiPRSV(f,T,P,stp)        
            print(Psat,"Psat found for steam",T)
            break
        P=P+stp
    return Psat


######################################################################
    
####fndpsatwith threads
def fndpsat_threadf(Pi,T,f,stp,ne):
    global counter_zeros, zero
    counter_zeros=0
    zero=0.0
    err=[0.0]*ne    
    for i in range (ne):
       
        if counter_zeros==1:
            break
        yi=Reynoldscolonna(Pi,T,f)
        err[i]=yi        
        
        if i>0:
            if err[i]*err[i-1]<0:
                gss=Pi-0.5*stp
                zero=fsolve (Reynoldscolonna,gss,(T,f))
                e=Reynoldscolonna(zero,T,f)
                if abs(e)>1.7e-06:
                    print("Psat not found",e)
                counter_zeros+=1
                break
        Pi+=stp 


def fnd_PsatThreat(f,T,Pi):   
           
    #    zero=0
    #    counter_zeros=0 
           
    nts=2#number of threads
    Pf=617#for water
    ne=1000
    net=int(ne/nts)+2
    stp=(Pf-Pi)/ne
    l=(Pf-Pi)
    #stp=stp/150
    
    threads=[]
    
    for i in range (nts):
        if i ==0:
            Pi=Pi
        if i>0:
            Pi+=l/nts
            
        thrd=threading.Thread(target=fndpsat_threadf,args=(Pi,T,f,stp,net))
        thrd.start()
        threads.append(thrd)
        
    for thrd in threads:
        thrd.join()
    
    Psat=zero
    return Psat

def fnd_psatiPRSV_thread(f,T):

    P=0.1*CP.PropsSI('Pcrit',f)
    #stp=-0.6e+05
    stp=-200000
    for i in range(10000000):
        msg=iPRSV_mas_PT(f,P,T)
        
        msg=msg[7]
        if msg=='Warning Two-phase zone':
            Pi=P
            Psat=fnd_PsatThreat(f,T,Pi)
            break
        P+=stp
    return Psat
###########################################################################    
###########################################################################
    
def fndpsat_Procc_f(Pi,T,f,stp,ne,zero):
    counter_zeros=0
    
    err=[0.0]*ne    
    for i in range (ne):
       
        if counter_zeros==1:
            break
        yi=Reynoldscolonna(Pi,T,f)
        err[i]=yi        
        
        if i>0:
            if err[i]*err[i-1]<0:
                gss=Pi-0.5*stp
                zero[0]=fsolve (Reynoldscolonna,gss,(T,f))
                e=Reynoldscolonna(zero[0],T,f)
                if abs(e)>1.7e-06:
                    print("Psat not found",e)
                break
        Pi+=stp 


def fnd_Psat_Process(f,T,Pi):
    nps=8 #NUMBER OF PROCESSES EMPLOYED
    zero=[0.0]*nps
    
    
    Pf=617#for water
    ne=1000
    ne_p=int(ne/nps)+2
    stp=(Pf-Pi)/ne
    l=(Pf-Pi)
    
    
    processes=[]
    for i in range (nps):
        if i ==0:
            Pi=Pi
        if i>0:
            Pi+=l/nps
        zero[i]=multiprocessing.Array('d',1)   
        process=multiprocessing.Process(target=fndpsat_Procc_f,args=(Pi,T,f,stp,ne_p,zero[i]))
        process.start()
        processes.append(process)
    
    
    for process in processes:
            process.join()
    
    res_zeros=[0.0]*nps
    for i in range (nps):
            res_zeros[i]=(zero[i][0])
            res_zeros[i]=float(res_zeros[i])
            if res_zeros[i]>0:
                Psat=res_zeros[i]
                break
    return Psat

def fnd_psatiPRSV_Process(f,T):
    P=0.2*CP.PropsSI('Pcrit',f)
    #stp=-0.6e+05
    stp=-20000
    for i in range(100000000):
        msg=iPRSV_mas_PT(f,P,T)
        msg=msg[7]
        if msg=='Warning Two-phase zone':
            Pi=P
            Psat=fnd_Psat_Process(f,T,Pi)
            e=Reynoldscolonna(Psat,T,f)
            if abs(e)>1.7e-07:
                    print("Psat not found",e)
            
            break
        P=P+stp
    return Psat


###################################################################
#################################################################






























