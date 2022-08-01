#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 18:29:10 2020

@author: andres
"""
from scipy.integrate import quad
import pandas as pd
import CoolProp.CoolProp as CP
import math
import numpy as np
from scipy.optimize import fsolve
from LibiPRSV import iPRSV_mas_PT, iPRSV_mas_PT_sat
from LibiPRSVsat import fnd_psatiPRSV, Ncl_rate, fnd_Tsat,fnd_psatiPRSV_thread,fnd_psatiPRSV_Process
from scipy.integrate import quad
from GeOsiPRSV import Moore_nozzle, Moses_nozzle
import matplotlib.pyplot as plt

def Ix(x,P4,P3,A4,A3):
        m=(P4-P3)/(A4-A3)
        b=P3-m*A3
        y=m*x+b       
        return y

def sft(T,f):
    Tc=CP.PropsSI('Tcrit',f)
    B=235.8/1000
    b=-0.625
    miu=1.256
    sft=B*((Tc-T)/Tc)**miu*(1+b*(Tc-T)/Tc)
    return sft







############################################################################################

def Chinos_M(xgss,P1,rho1,h1,u1,rcrit1,N,r,xi,xf,f,i,fs):
      c=i+1
       
      Nsum=0
      Bsum=0
      for j in range(c):
            if j<(c-1):
                Bsum+=r[j]*N[j]
                Nsum+=N[j]
            if j==(c-1):
                Bsum+=rcrit1*N[j]
                Nsum+=N[j]
    
      B=Bsum/Nsum  
      A=(xf-xi)/u1/fs
      
      P2=xgss[0]
      T2=xgss[1]
      
      "Gas properties"
      propg=iPRSV_mas_PT(f,P2,T2)
      P=propg[0]
      rhog=1/propg[1]
      Tg=propg[2]
      hg=propg[3]
      sg=propg[4]
      gamma=propg[8]
      
      cpg=propg[10]
      
      
      "Liquid properties"
      Tsat=fnd_Tsat(P,f)
      propsat_p=iPRSV_mas_PT_sat(f,P,Tsat)
      rhol=1/propsat_p[6]
      hl=propsat_p[7]
      sl=propsat_p[8]
      
      L=hg-hl
      delT=Tsat-Tg

      "Dropleth growth"
      R=8.31462/CP.PropsSI('molar_mass',f)
      
      sf=sft(Tg,f)
      C=P/(L*rhol*math.sqrt(2*math.pi*R*Tg))*(gamma+1)/(2*gamma)*cpg
      
      coeff=np.zeros(3)
      coeff[0]=1/(A*C)
      coeff[1]=-(B/(A*C)+delT)
      coeff[2]=(2*sf*Tsat)/(rhol*L)
      rm=np.roots(coeff)
 
      rm=max(rm)
      rm=rm.real
      dr_dt=C*(delT-(2*sf*Tsat)/(rhol*L*rm))
      delr=dr_dt*A
      delr=delr.real
      
      rnew=np.zeros(c)
      for j in range(c):
            if j<(c-1):
                rnew[j]=r[j]+delr
            if j==(c-1):
                rnew[j]=rcrit1+delr
      sumNrrr=0
      for j in range(c):
           sumNrrr+=N[j]*rnew[j]**(3.0)
      
      Y=4/3*math.pi*rhol*sumNrrr
      X=1-Y
        
      rho2=rhog/X
      h2=hg*X+hl*Y
    
        
      e1=h1+0.5*u1**2.0
        
      u2=math.sqrt(2*(e1-h2))
        
      A1=Moore_nozzle(xi)
      mass1=rho1*A1*u1
      A2=Moore_nozzle(xf)
      mass2=rho2*A2*u2
      
      I=quad(Ix,A1,A2,(P,P1,A2,A1))
      I=I[0]
      
      mom1=(P1+rho1*u1**2.0)*A1+I
      mom2=(P+rho2*u2**2.0)*A2
      ymass=mass2-mass1
      ymom=mom2-mom1
      diff=np.zeros(2)
      diff[0]=ymass
      diff[1]=ymom
      
      return diff

def Chinos_res(xgss,P1,rho1,h1,u1,rcrit1,N,r,xi,xf,f,i,fs):
      c=i+1
       
      Nsum=0
      Bsum=0
      for j in range(c):
            if j<(c-1):
                Bsum+=r[j]*N[j]
                Nsum+=N[j]
            if j==(c-1):
                Bsum+=rcrit1*N[j]
                Nsum+=N[j]
    
      B=Bsum/Nsum  
      A=(xf-xi)/u1/fs
       
      P2=xgss[0]
      T2=xgss[1]
      
      "Gas properties"
      propg=iPRSV_mas_PT(f,P2,T2)
      P=propg[0]
      rhog=1/propg[1]
      Tg=propg[2]
      hg=propg[3]
      sg=propg[4]
      gamma=propg[8]
      cpg=propg[10]
      
      
      "Liquid properties"
      Tsat=fnd_Tsat(P,f)
      propsat_p=iPRSV_mas_PT_sat(f,P,Tsat)
      rhol=1/propsat_p[6]
      hl=propsat_p[7]
      sl=propsat_p[8]
      
      L=hg-hl
      delT=Tsat-Tg
      
      "Dropleth growth"
      R=8.31462/CP.PropsSI('molar_mass',f)
      sf=sft(Tg,f)
      C=P/(L*rhol*math.sqrt(2*math.pi*R*Tg))*(gamma+1)/(2*gamma)*cpg
      coeff=np.zeros(3)
      coeff[0]=1/(A*C)
      coeff[1]=-(B/(A*C)+delT)
      coeff[2]=(2*sf*Tsat)/(rhol*L)
      rm=np.roots(coeff)
      rm=max(rm).real
      dr_dt=C*(delT-(2*sf*Tsat)/(rhol*L*rm))
      delr=dr_dt*A
      delr=delr.real
      #print(L/1000,": L")
      rnew=np.zeros(c)
      for j in range(c):
            if j<(c-1):
                rnew[j]=r[j]+delr
            if j==(c-1):
                rnew[j]=rcrit1+delr
      sumNrrr=0
      for j in range(c):
           sumNrrr+=N[j]*rnew[j]**(3.0)
      
      Y=4/3*math.pi*rhol*sumNrrr
      X=1-Y
        
      rho2=rhog/X
      h2=hg*X+hl*Y
      s2=sg*X+sl*Y
        
      e1=h1+0.5*u1**2.0
        
      u2=math.sqrt(2*(e1-h2))
        
      A1=Moore_nozzle(xi)
      mass1=rho1*A1*u1
      A2=Moore_nozzle(xf)
      mass2=rho2*A2*u2
      
      I=quad(Ix,A1,A2,(P,P1,A2,A1))
      I=I[0]
      mom1=(P1+rho1*u1**2.0)*A1+I
      mom2=(P+rho2*u2**2.0)*A2
      ymass=mass2-mass1
      ymom=mom2-mom1
      
      diff=np.zeros(2)
      diff[0]=ymass
      diff[1]=ymom
      
      return P,rho2,Tg,h2,s2,u2,Y,delr,rnew,rm
  
def Chinos_growth(P3,T3,u3,A2,grid,xcond,f,k):
    fs=10
    prop3=iPRSV_mas_PT(f,P3,T3)
    rho3=1/prop3[1]
    h3=prop3[3]
    s3=prop3[4]
    
    T3sat=fnd_Tsat(P3,f)
    P3sat=fnd_psatiPRSV_Process(f,T3)
    print(P3/P3sat)
    J3=Ncl_rate(f,P3,T3,P3sat,T3sat,k)
    rcrit3=J3[1]
    J3=J3[0]
    
    xf=7
    stp=(xf-xcond)/grid
    xi=xcond
    xf=xi+stp
    
    N3=J3*(xf-xi)/(rho3*u3)/fs
    
    Pi=np.zeros(grid)
    rhoi=np.zeros(grid)
    Ti=np.zeros(grid)
    hi=np.zeros(grid)
    si=np.zeros(grid)
    Si=np.zeros(grid)
    ui=np.zeros(grid)
    Yi=np.zeros(grid)
    SubCi=np.zeros(grid)
    Tsati=np.zeros(grid)
    Psati=np.zeros(grid)
    xii=np.zeros(grid)
    xif=np.zeros(grid)
    
    
    Ji=np.zeros(grid)
    Ni=np.zeros(grid)
    rcriti=np.zeros(grid)
    Jmax_counter=0
    ri=0
    rmi=np.zeros(grid)
    rnew=np.zeros(grid)
    delri=np.zeros(grid)
    rgroups=[0.0]*grid
    rgroup1=np.zeros(grid)
    rgroup2=np.zeros(grid)
    em=np.zeros(grid)
    emom=np.zeros(grid)
    ne=grid-1
    xgss=np.zeros(2)
    
    for i in range (ne):
        if i==0:
            Pi[i]=P3
            rhoi[i]=rho3
            Ti[i]=T3
            hi[i]=h3
            si[i]=s3
            ui[i]=u3
            Tsati[i]=T3sat
            xii[i]=xi
            xif[i]=xf
            SubCi[i]=Tsati[i]-Ti[i]
            Ji[i]=J3
            Ni[i]=N3
            rcriti[i]=rcrit3
            Psati[i]=P3sat
            Si[i]=Pi[i]/Psati[i]
            
        if xii[i]<2:
            xgss[0]=1.01*Pi[i]
            xgss[1]=1.01*Ti[i]
        else:        
           xgss[0]=0.9*Pi[i]
           xgss[1]=0.9*Ti[i]
            
            

        sol=fsolve(Chinos_M,xgss,(Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,i,fs))
        #sol=Chinos_M(xgss,Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,i,fs)
        diff=Chinos_M(sol,Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,i,fs)
           
        res=Chinos_res(sol,Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,i,fs)
        em[i+1]=diff[0]
        emom[i+1]=diff[1]
        Pi[i+1]=res[0]
        Ti[i+1]=res[2]
        rhoi[i+1]=res[1]
        Tsati[i+1]=fnd_Tsat(Pi[i+1],f)
        Psati[i+1]=fnd_psatiPRSV_Process(f,Ti[i+1])
        SubCi[i+1]=Tsati[i+1]-Ti[i+1]
        Si[i+1]=Pi[i+1]/Psati[i+1]
        
        hi[i+1]=res[3]
        si[i+1]=res[4]
        ui[i+1]=res[5]
        Yi[i+1]=res[6]
        delri[i+1]=res[7]
        rmi[i+1]=res[9]
        ri=res[8]
        rgroups[i+1]=ri
        
        nrate=Ncl_rate(f,Pi[i+1],Ti[i+1],Psati[i+1],Tsati[i+1],k)
        Ji[i+1]=nrate[0]
        rcriti[i+1]=nrate[1]
        
        print(xii[i])
        
        if Jmax_counter==0:
            if Ji[i+1]<Ji[i]:
                Jmax_counter=1
                xc=xii[i]
        if Jmax_counter==1:
            rnew[i+1]=rmi[i+1]
        
        xii[i+1]=xif[i]
        xif[i+1]=xii[i+1]+stp
        Ni[i+1]=Ji[i+1]*(xif[i+1]-xii[i+1])/(rhoi[i+1]*ui[i+1])/fs
        
        if Ti[i+1]<=276:
            break
        
        if i>=20:
            rgroup1[i+1]=rgroups[i+1][20]
        if i>=50:
            rgroup2[i+1]=rgroups[i+1][50]
    
    plt.xlabel("Nozzle Axis")  
    plt.ylabel("Pressure (bar)") 
    #plt.xlim(xii[0],xii[i-1])
    #plt.ylim(5000,13000)
    plt.plot(xii,Pi,"--")
    plt.grid()
      
    return ne,xii,Pi,Ti,ui,Yi,si,Si,Ji,rmi,rnew,em,emom,xc,rgroup1,rgroup2,SubCi,rcriti




def Gyarmathy_M_Moore(xgss,P1,rho1,h1,u1,rcrit1,N,r,xi,xf,f,i,fs):
      c=i+1
       
      Nsum=0
      Dsum=0
      for j in range(c):
            if j<(c-1):
                Dsum+=r[j]*N[j]
                Nsum+=N[j]
            if j==(c-1):
                Dsum+=rcrit1*N[j]
                Nsum+=N[j]
    
      D=Dsum/Nsum  
      C=(xf-xi)/u1/fs
      
      P2=xgss[0]
      T2=xgss[1]
      
      "Gas properties"
      propg=iPRSV_mas_PT(f,P2,T2)
      P=propg[0]
      rhog=1/propg[1]
      Tg=propg[2]
      hg=propg[3]
      sg=propg[4]
      gamma=propg[8]
      cpg=propg[10]
      
      mug=CP.PropsSI("viscosity","T",T2,"P",P2,f)
      lamda=CP.PropsSI("conductivity","T",T2,"P",P2,f)
      
      "Liquid properties"
      Tsat=fnd_Tsat(P,f)
      propsat_p=iPRSV_mas_PT_sat(f,P,Tsat)
      rhol=1/propsat_p[6]
      hl=propsat_p[7]
      sl=propsat_p[8]
      
      L=hg-hl
      delT=Tsat-Tg

      "Dropleth growth"
      
      R=8.31462/CP.PropsSI('molar_mass',f)
      A=1.5*mug*math.sqrt(R*Tg)/(2*P)
      B=lamda*delT/rhol/L
      coeff=np.zeros(4)
    
      coeff[0]=1
      coeff[1]=3.18*A-D
      coeff[2]=-(3.18*A*D+B*C)
      coeff[3]=B*C*rcrit1
    
    
      rm=np.roots(coeff)
      rm=max(rm)
      rm=rm.real
    
      drdt=B/(1+3.18*A/rm)*(rm-rcrit1)/rm**2.0    
      delr=rm-D  
      delr=delr.real
 

      
      rnew=np.zeros(c)
      for j in range(c):
            if j<(c-1):
                rnew[j]=r[j]+delr
            if j==(c-1):
                rnew[j]=rcrit1+delr
      sumNrrr=0
      for j in range(c):
           sumNrrr+=N[j]*rnew[j]**(3.0)
      
      Y=4/3*math.pi*rhol*sumNrrr
      X=1-Y
        
      rho2=rhog/X
      h2=hg*X+hl*Y
    
        
      e1=h1+0.5*u1**2.0
        
      u2=math.sqrt(2*(e1-h2))
        
      A1=Moore_nozzle(xi)
      mass1=rho1*A1*u1
      A2=Moore_nozzle(xf)
      mass2=rho2*A2*u2
      
      I=quad(Ix,A1,A2,(P,P1,A2,A1))
      I=I[0]
      
      mom1=(P1+rho1*u1**2.0)*A1+I
      mom2=(P+rho2*u2**2.0)*A2
      ymass=mass2-mass1
      ymom=mom2-mom1
      diff=np.zeros(2)
      diff[0]=ymass
      diff[1]=ymom
      
      return diff

def Gyarmathy_res_Moore(xgss,P1,rho1,h1,u1,rcrit1,N,r,xi,xf,f,i,fs):
      c=i+1
       
      Nsum=0
      Dsum=0
      for j in range(c):
            if j<(c-1):
                Dsum+=r[j]*N[j]
                Nsum+=N[j]
            if j==(c-1):
                Dsum+=rcrit1*N[j]
                Nsum+=N[j]
    
      D=Dsum/Nsum  
      C=(xf-xi)/u1/fs
      
      P2=xgss[0]
      T2=xgss[1]
      
      "Gas properties"
      propg=iPRSV_mas_PT(f,P2,T2)
      P=propg[0]
      rhog=1/propg[1]
      Tg=propg[2]
      hg=propg[3]
      sg=propg[4]
      cg=propg[6]
      gamma=propg[8]
      cpg=propg[10]
      
      mug=CP.PropsSI("viscosity","T",T2,"P",P2,f)
      lamda=CP.PropsSI("conductivity","T",T2,"P",P2,f)
      
      "Liquid properties"
      Tsat=fnd_Tsat(P,f)
      propsat_p=iPRSV_mas_PT_sat(f,P,Tsat)
      rhol=1/propsat_p[6]
      hl=propsat_p[7]
      sl=propsat_p[8]
      cl=propsat_p[11]
      
      print(cl)
      L=hg-hl
      delT=Tsat-Tg

      "Dropleth growth"
      
      R=8.31462/CP.PropsSI('molar_mass',f)
      A=1.5*mug*math.sqrt(R*Tg)/(2*P)
      B=lamda*delT/rhol/L
      coeff=np.zeros(4)
    
      coeff[0]=1
      coeff[1]=3.18*A-D
      coeff[2]=-(3.18*A*D+B*C)
      coeff[3]=B*C*rcrit1
    
    
      rm=np.roots(coeff)
      rm=max(rm)
      rm=rm.real
    
      drdt=B/(1+3.18*A/rm)*(rm-rcrit1)/rm**2.0    
      delr=rm-D  
      delr=delr.real
      Kn=1.5*mug*math.sqrt(R*Tg)/(2*rm*P)
 

      
      rnew=np.zeros(c)
      for j in range(c):
            if j<(c-1):
                rnew[j]=r[j]+delr
            if j==(c-1):
                rnew[j]=rcrit1+delr
      sumNrrr=0
      for j in range(c):
           sumNrrr+=N[j]*rnew[j]**(3.0)
      
      Y=4/3*math.pi*rhol*sumNrrr
      X=1-Y
        
      rho2=rhog/X
      h2=hg*X+hl*Y
      s2=sg*X+sl*Y
      c2=cg*X+cl*Y
    
        
      e1=h1+0.5*u1**2.0
        
      u2=math.sqrt(2*(e1-h2))
        
      A1=Moore_nozzle(xi)
      mass1=rho1*A1*u1
      A2=Moore_nozzle(xf)
      mass2=rho2*A2*u2
      
      I=quad(Ix,A1,A2,(P,P1,A2,A1))
      I=I[0]
      
      mom1=(P1+rho1*u1**2.0)*A1+I
      mom2=(P+rho2*u2**2.0)*A2
      ymass=mass2-mass1
      ymom=mom2-mom1
      diff=np.zeros(2)
      diff[0]=ymass
      diff[1]=ymom
      
      M2=u2/c2
      
      return P,rho2,Tg,h2,s2,u2,Y,delr,rnew,rm,Nsum,Kn,M2 




def Gyarmathy_growth_Moore(P3,T3,u3,A2,grid,xcond,f,k):
    fs=10
    prop3=iPRSV_mas_PT(f,P3,T3)
    rho3=1/prop3[1]
    h3=prop3[3]
    s3=prop3[4]
    M3=u3/prop3[6]
    
    T3sat=fnd_Tsat(P3,f)
    P3sat=fnd_psatiPRSV_Process(f,T3)
    print(P3/P3sat)
    J3=Ncl_rate(f,P3,T3,P3sat,T3sat,k)
    rcrit3=J3[1]
    J3=J3[0]
    
    xf=7
    stp=(xf-xcond)/grid
    xi=xcond
    xf=xi+stp
    
    N3=J3*(xf-xi)/(rho3*u3)/fs
    
    Pi=np.zeros(grid)
    rhoi=np.zeros(grid)
    Ti=np.zeros(grid)
    hi=np.zeros(grid)
    si=np.zeros(grid)
    Si=np.zeros(grid)
    ui=np.zeros(grid)
    Yi=np.zeros(grid)
    SubCi=np.zeros(grid)
    Tsati=np.zeros(grid)
    Psati=np.zeros(grid)
    Mi=np.zeros(grid)
    eener=np.zeros(grid)
    
    xii=np.zeros(grid)
    xif=np.zeros(grid)
    
    
    Ji=np.zeros(grid)
    Ni=np.zeros(grid)
    Nsum=np.zeros(grid)
    rcriti=np.zeros(grid)
    Jmax_counter=0
    ri=0
    rmi=np.zeros(grid)
    rnew=np.zeros(grid)
    delri=np.zeros(grid)
    ei=np.zeros(grid)
    Kn=np.zeros(grid)
    
    rgroups=[0.0]*grid
    rgroup1=np.zeros(grid)
    rgroup2=np.zeros(grid)
    em=np.zeros(grid)
    emom=np.zeros(grid)
    ne=grid-1
    xgss=np.zeros(2)
    
    for i in range (ne):
        if i==0:
            Pi[i]=P3
            rhoi[i]=rho3
            Ti[i]=T3
            hi[i]=h3
            si[i]=s3
            ui[i]=u3
            Tsati[i]=T3sat
            xii[i]=xi
            xif[i]=xf
            SubCi[i]=Tsati[i]-Ti[i]
            Ji[i]=J3
            Ni[i]=N3
            Nsum[i]=N3
            Mi[i]=M3
            rcriti[i]=rcrit3
            Psati[i]=P3sat
            ei[i]=h3+0.5*u3**2.0
            Si[i]=Pi[i]/Psati[i]
            
        if xii[i]<2:
            xgss[0]=1.001*Pi[i]
            xgss[1]=1.001*Ti[i]
        else:        
           xgss[0]=0.98*Pi[i]
           xgss[1]=0.98*Ti[i]
            
            

        sol=fsolve(Gyarmathy_M_Moore,xgss,(Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,i,fs),xtol=1.5e-20)
        #sol=Chinos_M(xgss,Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,i,fs)
        diff=Gyarmathy_M_Moore(sol,Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,i,fs)
           
        res=Gyarmathy_res_Moore(sol,Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,i,fs)
        em[i+1]=diff[0]
        emom[i+1]=diff[1]
        Pi[i+1]=res[0]
        Ti[i+1]=res[2]
        rhoi[i+1]=res[1]
        Tsati[i+1]=fnd_Tsat(Pi[i+1],f)
        Psati[i+1]=fnd_psatiPRSV_Process(f,Ti[i+1])
        SubCi[i+1]=Tsati[i+1]-Ti[i+1]
        Si[i+1]=Pi[i+1]/Psati[i+1]
        emom[i+1]=diff[1]
        em[i+1]=diff[0]
        
        hi[i+1]=res[3]
        si[i+1]=res[4]
        ui[i+1]=res[5]
        Yi[i+1]=res[6]
        delri[i+1]=res[7]
        rmi[i+1]=res[9]
        Mi[i+1]=res[12]
        Nsum[i+1]=res[10]
        ri=res[8]
        rgroups[i+1]=ri
        ei[i+1]=hi[i+1]+0.5*ui[i+1]**2.0
        eener[i+1]=ei[i+1]-ei[i]
        nrate=Ncl_rate(f,Pi[i+1],Ti[i+1],Psati[i+1],Tsati[i+1],k)
        Ji[i+1]=nrate[0]
        rcriti[i+1]=nrate[1]
        Kn[i+1]=res[11]
        
        print(xii[i],Ti[i],Ji[i])
        
        if Jmax_counter==0:
            if Ji[i+1]<Ji[i]:
                Jmax_counter=1
                xc=xii[i]
                Pw=Pi[i]
                sw=si[i]
                Yw=Yi[i]
                Tw=Ti[i]
        if Jmax_counter==1:
            rnew[i+1]=rmi[i+1]
        
        xii[i+1]=xif[i]
        xif[i+1]=xii[i+1]+stp
        Ni[i+1]=Ji[i+1]*(xif[i+1]-xii[i+1])/(rhoi[i+1]*ui[i+1])/fs
        
        if Ti[i+1]<=276:
            break
        
        if i>=20:
            rgroup1[i+1]=rgroups[i+1][20]
        if i>=50:
            rgroup2[i+1]=rgroups[i+1][50]
    
    plt.xlabel("Nozzle Axis")  
    plt.ylabel("Pressure (bar)") 
    #plt.xlim(xii[0],xii[i-1])
    #plt.ylim(5000,13000)
    plt.plot(xii,Pi,"--")
    plt.grid()
    
    
    return ne,xii,Pi,Ti,ui,Yi,si,Si,Ji,rmi,rnew,em,emom,xc,rgroup1,rgroup2,SubCi,rcriti,delri,Nsum,Kn,Pw,Tw,Yw,eener,ei,Mi 

##################################################################################################
############################################################################


def Gyarmathy_M_Moses(xgss,P1,rho1,h1,u1,rcrit1,N,r,xi,xf,f,i,fs):
      c=i+1
       
      Nsum=0
      Dsum=0
      for j in range(c):
            if j<(c-1):
                Dsum+=r[j]*N[j]
                Nsum+=N[j]
            if j==(c-1):
                Dsum+=rcrit1*N[j]
                Nsum+=N[j]
    
      D=Dsum/Nsum  
      C=(xf-xi)/u1/fs
      
      P2=xgss[0]
      T2=xgss[1]
      
      "Gas properties"
      propg=iPRSV_mas_PT(f,P2,T2)
      P=propg[0]
      rhog=1/propg[1]
      Tg=propg[2]
      hg=propg[3]
      sg=propg[4]
      gamma=propg[8]
      cpg=propg[10]
      
      mug=CP.PropsSI("viscosity","T",T2,"P",P2,f)
      lamda=CP.PropsSI("conductivity","T",T2,"P",P2,f)
      
      "Liquid properties"
      Tsat=fnd_Tsat(P,f)
      propsat_p=iPRSV_mas_PT_sat(f,P,Tsat)
      rhol=1/propsat_p[6]
      hl=propsat_p[7]
      sl=propsat_p[8]
      
      
      L=hg-hl
      delT=Tsat-Tg

      "Dropleth growth"
      
      R=8.31462/CP.PropsSI('molar_mass',f)
      A=1.5*mug*math.sqrt(R*Tg)/(2*P)
      B=lamda*delT/rhol/L
      coeff=np.zeros(4)
    
      coeff[0]=1
      coeff[1]=3.18*A-D
      coeff[2]=-(3.18*A*D+B*C)
      coeff[3]=B*C*rcrit1
    
    
      rm=np.roots(coeff)
      rm=max(rm)
      rm=rm.real
    
      drdt=B/(1+3.18*A/rm)*(rm-rcrit1)/rm**2.0    
      delr=rm-D  
      delr=delr.real
 

      
      rnew=np.zeros(c)
      for j in range(c):
            if j<(c-1):
                rnew[j]=r[j]+delr
            if j==(c-1):
                rnew[j]=rcrit1+delr
      sumNrrr=0
      for j in range(c):
           sumNrrr+=N[j]*rnew[j]**(3.0)
      
      Y=4/3*math.pi*rhol*sumNrrr
      X=1-Y
        
      rho2=rhog/X
      h2=hg*X+hl*Y
    
        
      e1=h1+0.5*u1**2.0
        
      u2=math.sqrt(2*(e1-h2))
        
      A1=Moses_nozzle(xi)
      mass1=rho1*A1*u1
      A2=Moses_nozzle(xf)
      mass2=rho2*A2*u2
      
      I=quad(Ix,A1,A2,(P,P1,A2,A1))
      I=I[0]
      
      mom1=(P1+rho1*u1**2.0)*A1+I
      mom2=(P+rho2*u2**2.0)*A2
      ymass=mass2-mass1
      ymom=mom2-mom1
      diff=np.zeros(2)
      diff[0]=ymass
      diff[1]=ymom
      
      return diff


def Gyarmathy_res_Moses(xgss,P1,rho1,h1,u1,rcrit1,N,r,xi,xf,f,i,fs):
      c=i+1
       
      Nsum=0
      Dsum=0
      for j in range(c):
            if j<(c-1):
                Dsum+=r[j]*N[j]
                Nsum+=N[j]
            if j==(c-1):
                Dsum+=rcrit1*N[j]
                Nsum+=N[j]
    
      D=Dsum/Nsum  
      C=(xf-xi)/u1/fs
      
      P2=xgss[0]
      T2=xgss[1]
      
      "Gas properties"
      propg=iPRSV_mas_PT(f,P2,T2)
      P=propg[0]
      rhog=1/propg[1]
      Tg=propg[2]
      hg=propg[3]
      sg=propg[4]
      cg=propg[6]
      gamma=propg[8]
      cpg=propg[10]
      
      mug=CP.PropsSI("viscosity","T",T2,"P",P2,f)
      lamda=CP.PropsSI("conductivity","T",T2,"P",P2,f)
      
      "Liquid properties"
      Tsat=fnd_Tsat(P,f)
      propsat_p=iPRSV_mas_PT_sat(f,P,Tsat)
      rhol=1/propsat_p[6]
      hl=propsat_p[7]
      sl=propsat_p[8]
      cl=propsat_p[11]
      Tl=Tsat
      
      print(cl)
      L=hg-hl
      delT=Tsat-Tg

      "Dropleth growth"
      
      R=8.31462/CP.PropsSI('molar_mass',f)
      A=1.5*mug*math.sqrt(R*Tg)/(2*P)
      B=lamda*delT/rhol/L
      coeff=np.zeros(4)
    
      coeff[0]=1
      coeff[1]=3.18*A-D
      coeff[2]=-(3.18*A*D+B*C)
      coeff[3]=B*C*rcrit1
    
    
      rm=np.roots(coeff)
      rm=max(rm)
      rm=rm.real
    
      drdt=B/(1+3.18*A/rm)*(rm-rcrit1)/rm**2.0    
      delr=rm-D  
      delr=delr.real
      Kn=1.5*mug*math.sqrt(R*Tg)/(2*rm*P)
 

      
      rnew=np.zeros(c)
      for j in range(c):
            if j<(c-1):
                rnew[j]=r[j]+delr
            if j==(c-1):
                rnew[j]=rcrit1+delr
      sumNrrr=0
      for j in range(c):
           sumNrrr+=N[j]*rnew[j]**(3.0)
      
      Y=4/3*math.pi*rhol*sumNrrr
      X=1-Y
        
      rho2=rhog/X
      h2=hg*X+hl*Y
      s2=sg*X+sl*Y
      c2=cg*X+cl*Y
    
        
      e1=h1+0.5*u1**2.0
        
      u2=math.sqrt(2*(e1-h2))
        
      A1=Moses_nozzle(xi)
      mass1=rho1*A1*u1
      A2=Moses_nozzle(xf)
      mass2=rho2*A2*u2
      
      I=quad(Ix,A1,A2,(P,P1,A2,A1))
      I=I[0]
      
      mom1=(P1+rho1*u1**2.0)*A1+I
      mom2=(P+rho2*u2**2.0)*A2
      ymass=mass2-mass1
      ymom=mom2-mom1
      diff=np.zeros(2)
      diff[0]=ymass
      diff[1]=ymom
      
      M2=u2/c2
      
      return P,rho2,Tg,h2,s2,u2,Y,delr,rnew,rm,Nsum,Kn,M2,Tl 

def Gyarmathy_growth_Moses(P3,T3,u3,A2,grid,xcond,f,k):
    fs=100
    prop3=iPRSV_mas_PT(f,P3,T3)
    rho3=1/prop3[1]
    h3=prop3[3]
    s3=prop3[4]
    M3=u3/prop3[6]
    
    T3sat=fnd_Tsat(P3,f)
    P3sat=fnd_psatiPRSV_Process(f,T3)
    print(P3/P3sat)
    J3=Ncl_rate(f,P3,T3,P3sat,T3sat,k)
    rcrit3=J3[1]
    J3=J3[0]
    
    xf=16
    stp=(xf-xcond)/grid
    xi=xcond
    xf=xi+stp
    
    N3=J3*(xf-xi)/(rho3*u3)/fs
    
    Pi=np.zeros(grid)
    rhoi=np.zeros(grid)
    Ti=np.zeros(grid)
    hi=np.zeros(grid)
    si=np.zeros(grid)
    Si=np.zeros(grid)
    ui=np.zeros(grid)
    Yi=np.zeros(grid)
    SubCi=np.zeros(grid)
    Tsati=np.zeros(grid)
    Psati=np.zeros(grid)
    Mi=np.zeros(grid)
    Tl=np.zeros(grid)
    Tlnew=np.zeros(grid)
    eener=np.zeros(grid)
    
    xii=np.zeros(grid)
    xif=np.zeros(grid)
    
    
    Ji=np.zeros(grid)
    Ni=np.zeros(grid)
    Nsum=np.zeros(grid)
    rcriti=np.zeros(grid)
    Jmax_counter=0
    ri=0
    rmi=np.zeros(grid)
    rnew=np.zeros(grid)
    delri=np.zeros(grid)
    ei=np.zeros(grid)
    Kn=np.zeros(grid)
    
    rgroups=[0.0]*grid
    rgroup1=np.zeros(grid)
    rgroup2=np.zeros(grid)
    em=np.zeros(grid)
    emom=np.zeros(grid)
    ne=grid-1
    xgss=np.zeros(2)
    
    for i in range (ne):
        if i==0:
            Pi[i]=P3
            rhoi[i]=rho3
            Ti[i]=T3
            hi[i]=h3
            si[i]=s3
            ui[i]=u3
            Tsati[i]=T3sat
            xii[i]=xi
            xif[i]=xf
            SubCi[i]=Tsati[i]-Ti[i]
            Ji[i]=J3
            Ni[i]=N3
            Nsum[i]=N3
            Mi[i]=M3
            rcriti[i]=rcrit3
            Psati[i]=P3sat
            ei[i]=h3+0.5*u3**2.0
            Si[i]=Pi[i]/Psati[i]
            
        if xii[i]<8.22:
            xgss[0]=1.001*Pi[i]
            xgss[1]=1.001*Ti[i]
        else:        
           xgss[0]=0.98*Pi[i]
           xgss[1]=0.98*Ti[i]
            
            

        sol=fsolve(Gyarmathy_M_Moses,xgss,(Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,i,fs),xtol=1.5e-15)
        #sol=Chinos_M(xgss,Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,i,fs)
        diff=Gyarmathy_M_Moses(sol,Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,i,fs)
           
        res=Gyarmathy_res_Moses(sol,Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,i,fs)
        em[i+1]=diff[0]
        emom[i+1]=diff[1]
        Pi[i+1]=res[0]
        Ti[i+1]=res[2]
        rhoi[i+1]=res[1]
        Tsati[i+1]=fnd_Tsat(Pi[i+1],f)
        Psati[i+1]=fnd_psatiPRSV_Process(f,Ti[i+1])
        SubCi[i+1]=Tsati[i+1]-Ti[i+1]
        Si[i+1]=Pi[i+1]/Psati[i+1]
        emom[i+1]=diff[1]
        em[i+1]=diff[0]
        
        hi[i+1]=res[3]
        si[i+1]=res[4]
        ui[i+1]=res[5]
        Yi[i+1]=res[6]
        delri[i+1]=res[7]
        rmi[i+1]=res[9]
        Mi[i+1]=res[12]
        Tl[i+1]=res[13]
        Nsum[i+1]=res[10]
        ri=res[8]
        rgroups[i+1]=ri
        ei[i+1]=hi[i+1]+0.5*ui[i+1]**2.0
        eener[i+1]=ei[i+1]-ei[i]
        nrate=Ncl_rate(f,Pi[i+1],Ti[i+1],Psati[i+1],Tsati[i+1],k)
        Ji[i+1]=nrate[0]
        rcriti[i+1]=nrate[1]
        Kn[i+1]=res[11]
        
        print(xii[i],Ti[i],Ji[i])
        
        if Jmax_counter==0:
            if Ji[i+1]<Ji[i]:
                if(xii[i])>8.22: 
                
                    Jmax_counter=1
                    xc=xii[i]
                    Pw=Pi[i]
                    sw=si[i]
                    Yw=Yi[i]
                    Tw=Ti[i]
                    c=i
        if Jmax_counter==1:
            Tlnew[i+1]=Tl[i+1]
            rnew[i+1]=rmi[i+1]
            
            if i>=c-60:
              rgroup1[i+1]=rgroups[i+1][c-60]
            if i>=c+30*4:
              rgroup2[i+1]=rgroups[i+1][c+30*4]
        
        xii[i+1]=xif[i]
        xif[i+1]=xii[i+1]+stp
        Ni[i+1]=Ji[i+1]*(xif[i+1]-xii[i+1])/(rhoi[i+1]*ui[i+1])/fs
        
        if Ti[i+1]<=276:
            break
        """
        if i>=20:
            rgroup1[i+1]=rgroups[i+1][20]
        if i>=50:
            rgroup2[i+1]=rgroups[i+1][50]
        """
    """
    plt.xlabel("Nozzle Axis")  
    plt.ylabel("Pressure (bar)") 
    #plt.xlim(xii[0],xii[i-1])
    #plt.ylim(5000,13000)
    plt.plot(xii,Pi,"--")
    plt.grid()
    """
    
    return ne,xii,Pi,Ti,ui,Yi,si,Si,Ji,rmi,rnew,em,emom,xc,rgroup1,rgroup2,SubCi,rcriti,delri,Nsum,Kn,Pw,Tw,Yw,eener,ei,Mi,Tl,Tlnew 




##Moore##

"""
f="Water"
P3=16767.03327794258
T3=324.7144704796638
u3=354.2870350106649
xcond=1.3749999999999927
grid=200
k=0.5
A2=1



drpl=Gyarmathy_growth_Moore(P3,T3,u3,A2,grid,xcond,f,k)
"""
    









    
    