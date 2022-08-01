#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 15:47:40 2020

@author: andres
"""

from Geos import B4M
from Coollib import Cool_PT, Cool_Ps, Cool_hs, Cool_Psat
from Coollib import Cool_Tsat, Cool_PT, Cool_Psat, Cool_Tsat
from Cool_sat import T_sat, p_sat
import math
import CoolProp.CoolProp as CP
import numpy as np
from scipy.optimize import fsolve
from Cool_sat import Ncl_rate,Ncl_rate3 
import matplotlib.pyplot as plt
from scipy.integrate import quad
import pandas as pd

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



def Gyarmathy_M(xgss,P1,rho1,h1,u1,rcrit1,N,r,xi,xf,f,EoS,i):
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
    C=(xf-xi)/u1/1000
  
        
    
    P2=xgss[0]
    T2=xgss[1]
    
    "Gas properties"
    propg=Cool_PT(P2,T2,f,EoS)
    P=propg[0]
    rhog=propg[1]
    Tg=propg[2]
    hg=propg[3]
    sg=propg[4]

    mug=propg[11]
    
    "Liquid properties"
    Tsat=T_sat(P,f,EoS)
    propl=Cool_Tsat(Tsat,f,EoS)
    rhol=1/propl[6]
    hl=propl[7]
    sl=propl[8]
    
    lamda=propg[10]
    
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
    
    drdt=B/(1+3.18*A/rm)*(rm-rcrit1)/rm**2.0    
    delr=rm-D  
    
    
    rnew=np.zeros(c)    
    for j in range(c):
        if j<(c-1):
            rnew[j]=r[j]+delr
        if j==(c-1):
            rnew[j]=rcrit1+delr
    
    sumNr=0
    for j in range(c):
        sumNr+=N[j]*rnew[j]**(3.0)
        


    Y=4/3*math.pi*rhol*sumNr
    X=1-Y
    
    rho2=rhog/X
    h2=hg*X+hl*Y
    s2=sg*X+sl*Y
    
    e1=h1+0.5*u1**2.0
    
    u2=math.sqrt(2*(e1-h2))
    
    A1=B4M(xi)
    mass1=rho1*A1*u1
    A2=B4M(xf)
    mass2=rho2*A2*u2
    
    
    I=quad(Ix,A1,A2,(P,P1,A2,A1))
    I=I[0]
    
    mom1=(P1+rho1*u1**2.0)*A1+I
    mom2=(P+rho2*u2**2.0)*A2
    
    ymass=mass2-mass1
    ymom=mom2-mom1
    return ymass,ymom

def Gyarmathy_res(xgss,P1,rho1,h1,u1,rcrit1,N,r,xi,xf,f,EoS,i):
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
    C=(xf-xi)/u1/1000
  
    P2=xgss[0]
    T2=xgss[1]
    
    "Gas properties"
    propg=Cool_PT(P2,T2,f,EoS)
    P=propg[0]
    rhog=propg[1]
    Tg=propg[2]
    hg=propg[3]
    sg=propg[4]
    gamma=propg[7]
    cpg=propg[9]
    cg=propg[6]
    mug=propg[11]
    
    "Liquid properties"
    Tsat=T_sat(P,f,EoS)
    propl=Cool_Tsat(Tsat,f,EoS)
    rhol=1/propl[6]
    hl=propl[7]
    sl=propl[8]
    cl=propl[11]
    lamda=propg[10]
    Tl=Tsat
    
    
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
    
    drdt=B/(1+3.18*A/rm)*(rm-rcrit1)/rm**2.0    
    delr=rm-D  
    Kn=1.5*mug*math.sqrt(R*Tg)/(2*rm*P)
    
    rnew=np.zeros(c)    
    for j in range(c):
        if j<(c-1):
            rnew[j]=r[j]+delr
        if j==(c-1):
            rnew[j]=rcrit1+delr
    
    sumNr=0
    for j in range(c):
        sumNr+=N[j]*rnew[j]**(3.0)
        

        
    Y=4/3*math.pi*rhol*sumNr
    X=1-Y
    
    rho2=rhog/X
    h2=hg*X+hl*Y
    s2=sg*X+sl*Y
    c2=cg*X+cl*Y
    e1=h1+0.5*u1**2.0
    
    u2=math.sqrt(2*(e1-h2))
    A1=B4M(xi)
    mass1=rho1*A1*u1
    A2=B4M(xf)
    mass2=rho2*A2*u2
       
    I=quad(Ix,A1,A2,(P,P1,A2,A1))
    I=I[0]
    
    mom1=(P1+rho1*u1**2.0)*A1+I
    mom2=(P+rho2*u2**2.0)*A2
    M2=u2/c2
    

    return P2,rho2,Tg,h2,s2,u2,Y,delr,rnew,rm,Nsum,Kn,M2,Tl


def Gyarmathy_growth(P3,T3,u3,A2,grid,xcond,f,EoS,k):
    
    prop3=Cool_PT(P3,T3,f,EoS)
    s3=prop3[4]
    h3=prop3[3]
    M3=u3/prop3[6]
    rho3=prop3[1]   


    xf=30 
    stp=(xf-xcond)/grid
    xi=xcond
    xf=xi+stp
    
    J3=Ncl_rate(P3,T3,f,EoS,k)
    rcrit3=J3[1]
    J3=J3[0]
    N3=J3*(xf-xi)/(rho3*u3)/1000
    
    Pi=np.zeros(grid)
    rhoi=np.zeros(grid)
    Ti=np.zeros(grid)
    hi=np.zeros(grid)
    si=np.zeros(grid)
    Si=np.zeros(grid)
    ui=np.zeros(grid)
    Yi=np.zeros(grid)
    Mi=np.zeros(grid)
    SubCi=np.zeros(grid)
    
    xii=np.zeros(grid)
    xif=np.zeros(grid)
    
    xgss=np.zeros(2)
    Ji=np.zeros(grid)
    Ni=np.zeros(grid)
    Nsum=np.zeros(grid)
    ei=np.zeros(grid)
    Kn=np.zeros(grid)
    rcriti=np.zeros(grid)
    
    rmi=np.zeros(grid)
    rnew=np.zeros(grid)
    delri=np.zeros(grid)
    rgroups=[0.0]*grid
    rgroup1=np.zeros(grid)
    rgroup2=np.zeros(grid)
    em=np.zeros(grid)
    eener=np.zeros(grid)
    emom=np.zeros(grid)
    Tl=np.zeros(grid)
    Tlnew=np.zeros(grid)
    ne=grid-1
    
    Jmax_counter=0
    ri=0
    for i in range(ne):
        if i==0:
            
            Pi[i]=P3
            rhoi[i]=rho3
            Ti[i]=T3
            hi[i]=h3
            si[i]=s3
            ui[i]=u3
            ei[i]=h3+0.5*u3**2.0
            xii[i]=xi
            xif[i]=xf
            Mi[i]=M3
            Ji[i]=J3
            Ni[i]=N3
            Nsum[i]=N3
            rcriti[i]=rcrit3
            Si[i]=Pi[i]/(p_sat(Ti[i],f,EoS))
            SubCi[i]=T_sat(Pi[i],f,EoS)-Ti[i]
            xgss[0]=1.001*P3
            xgss[1]=1.01*T3
         
        xgss[0]=1.001*Pi[i]
        xgss[1]=1.01*Ti[i]  
       
            
            
        sol=fsolve(Gyarmathy_M,xgss,(Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,EoS,i),xtol=1.8e-015)
        diff=Gyarmathy_M(sol,Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,EoS,i)
        print(xii[i])
        res=Gyarmathy_res(sol,Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,EoS,i)
        em[i+1]=diff[0]
        emom[i+1]=diff[1]
        Pi[i+1]=res[0]
        rhoi[i+1]=res[1]
        Ti[i+1]=res[2]
        Si[i+1]=Pi[i+1]/(p_sat(Ti[i+1],f,EoS))
        SubCi[i+1]=T_sat(Pi[i+1],f,EoS)-Ti[i+1]
        Nsum[i+1]=res[10]
        hi[i+1]=res[3]
        si[i+1]=res[4]
        ui[i+1]=res[5]
        ei[i+1]=hi[i+1]+0.5*ui[i+1]**2.0
        eener[i+1]=ei[i+1]-ei[i]
        Yi[i+1]=res[6]
        Kn[i+1]=res[11]
        Mi[i+1]=res[12]
        delri[i+1]=res[7]
        rmi[i+1]=res[9]
        Tl[i+1]=res[13]
        
        ri=res[8]
        rgroups[i+1]=ri
        nrate=Ncl_rate(Pi[i+1],Ti[i+1],f,EoS,k)
        Ji[i+1]=nrate[0]
        
        if Jmax_counter==0:
            if Ji[i+1]<Ji[i]:
                print("shit",Pi[i])
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
            if i>=c+200:
              rgroup2[i+1]=rgroups[i+1][c+200]
              
        rcriti[i+1]=nrate[1]
        
        xii[i+1]=xif[i]
        xif[i+1]=xii[i+1]+stp
        Ni[i+1]=Ji[i+1]*(xif[i+1]-xii[i+1])/(rhoi[i+1]*ui[i+1])/1000
        
        

    """        
    plt.figure(1)       
    plt.xlabel("Nozzle Axis")  
    plt.ylabel("Pressure (bar)")      
    plt.plot(xii,Pi,"--")
    plt.grid()
    """
    """
    plt.figure(2)       
    plt.xlabel("Nozzle Axis")  
    plt.ylabel("Pressure (bar)")      
    plt.plot(xii,Ji,"--")
    plt.grid()
    
    plt.figure(3)       
    plt.xlabel("Nozzle Axis")  
    plt.ylabel("Pressure (bar)")      
    plt.plot(xii,rmi,"--")
    plt.grid()
    """

      
     
      
    return ne,xii,Pi,Ti,ui,Yi,si,Si,Ji,rmi,rnew,em,emom,xc,rgroup1,rgroup2,SubCi,rcriti,delri,Nsum,Kn,Pw,Tw,Yw,eener,ei,Mi,Tl,Tlnew


##########################################################
def Young_M(xgss,P1,rho1,h1,u1,rcrit1,N,r,xi,xf,f,EoS,i):
     c=i+1
     Nsum=0
     Esum=0
     for j in range(c):
         if j<(c-1):
             Esum+=r[j]*N[j]
             Nsum+=N[j]
         if j==(c-1):
             Esum+=rcrit1*N[j]
             Nsum+=N[j]
     E=Esum/Nsum
     F=(xf-xi)/u1/1000
     
     P2=xgss[0]
     T2=xgss[1]
     
     "Gas properties"
     propg=Cool_PT(P2,T2,f,EoS)
     P=propg[0]
     rhog=propg[1]
     Tg=propg[2]
     hg=propg[3]
     sg=propg[4]
     gamma=propg[7]
     cpg=propg[9]
     lamda=propg[10]
     mug=propg[11]
     Prg=cpg*mug/lamda
     
     propsat_p=Cool_Psat(P,f,EoS)
     Tsat=propsat_p[1]
     rhol=1/propsat_p[6]
     hl=propsat_p[7]
     sl=propsat_p[8]
     
     L=hg-hl
     R=8.31462/CP.PropsSI('molar_mass',f)
     
     "Droplet growth"
     alpha=0
     theta=R*Tsat/L*(alpha-0.5-0.5*(gamma+1)/(2*gamma)*(cpg*Tsat/L))
     delT=Tsat-Tg
     
     A=lamda*delT
     B=3.78*(1-theta)/Prg
     C=rhol*L
     D=1.5*mug*math.sqrt(R*Tg)/(2*P)
     
     coeff=np.zeros(4)
     coeff[0]=C
     coeff[1]=C*B*D-C*E
     coeff[2]=-(C*B*D*E+F*A)
     coeff[3]=F*A*rcrit1
     
     roots=np.roots(coeff)
     rm=max(roots)
     
     drdt=A/C*(1-rcrit1/rm)/(rm+B*D)
     delr=drdt*F
     
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
    
     A1=B4M(xi)
     mass1=rho1*A1*u1
     A2=B4M(xf)
     mass2=rho2*A2*u2
    
    
     I=quad(Ix,A1,A2,(P,P1,A2,A1))
     I=I[0]
    
     mom1=(P1+rho1*u1**2.0)*A1+I
     mom2=(P+rho2*u2**2.0)*A2
    
     ymass=mass2-mass1
     ymom=mom2-mom1
     return ymass,ymom
 
def Young_res(xgss,P1,rho1,h1,u1,rcrit1,N,r,xi,xf,f,EoS,i):
              
     c=i+1
     Nsum=0
     Esum=0
     for j in range(c):
         if j<(c-1):
             Esum+=r[j]*N[j]
             Nsum+=N[j]
         if j==(c-1):
             Esum+=rcrit1*N[j]
             Nsum+=N[j]
     E=Esum/Nsum
     F=(xf-xi)/u1/1000
     
     P2=xgss[0]
     T2=xgss[1]
     
     "Gas properties"
     propg=Cool_PT(P2,T2,f,EoS)
     P=propg[0]
     rhog=propg[1]
     Tg=propg[2]
     hg=propg[3]
     sg=propg[4]
     gamma=propg[7]
     cpg=propg[9]
     lamda=propg[10]
     mug=propg[11]
     Prg=cpg*mug/lamda
     
     propsat_p=Cool_Psat(P,f,EoS)
     Tsat=propsat_p[1]
     rhol=1/propsat_p[6]
     hl=propsat_p[7]
     sl=propsat_p[8]
     
     L=hg-hl
     R=8.31462/CP.PropsSI('molar_mass',f)
     
     "Droplet growth"
     alpha=0
     theta=R*Tsat/L*(alpha-0.5-0.5*(gamma+1)/(2*gamma)*(cpg*Tsat/L))
     delT=Tsat-Tg
     
     A=lamda*delT
     B=3.78*(1-theta)/Prg
     C=rhol*L
     D=1.5*mug*math.sqrt(R*Tg)/(2*P)
     
     coeff=np.zeros(4)
     coeff[0]=C
     coeff[1]=C*B*D-C*E
     coeff[2]=-(C*B*D*E+F*A)
     coeff[3]=F*A*rcrit1
     
     roots=np.roots(coeff)
     rm=max(roots)
     
     drdt=A/C*(1-rcrit1/rm)/(rm+B*D)
     delr=drdt*F
     
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
    
     A1=B4M(xi)
     mass1=rho1*A1*u1
     A2=B4M(xf)
     mass2=rho2*A2*u2
    
    
     I=quad(Ix,A1,A2,(P,P1,A2,A1))
     I=I[0]
    
     mom1=(P1+rho1*u1**2.0)*A1+I
     mom2=(P+rho2*u2**2.0)*A2
    
     ymass=mass2-mass1
     ymom=mom2-mom1
     return P2,rho2,Tg,h2,s2,u2,Y,delr,rnew,rm

    

def Young_growth(P3,T3,u3,A2,grid,xcond,f,EoS,k):

    prop3=Cool_PT(P3,T3,f,EoS)
    s3=prop3[4]
    h3=prop3[3]
    rho3=prop3[1]
    
    J3=Ncl_rate(P3,T3,f,EoS,k)
    rcrit3=J3[1]
    J3=J3[0]

    xf=30
    
   
    stp=(xf-xcond)/grid

    xi=xcond
    xf=xi+stp
    
    N3=J3*(xf-xi)/(rho3*u3)/1000
    
    Pi=np.zeros(grid)
    rhoi=np.zeros(grid)
    Ti=np.zeros(grid)
    hi=np.zeros(grid)
    si=np.zeros(grid)
    Si=np.zeros(grid)
    ui=np.zeros(grid)
    Yi=np.zeros(grid)
    
    
    xii=np.zeros(grid)
    xif=np.zeros(grid)
    
    xgss=np.zeros(2)
    Ji=np.zeros(grid)
    Ni=np.zeros(grid)
    rcriti=np.zeros(grid)
    ri=0
    rmi=np.zeros(grid)
    rnew=np.zeros(grid)
    delri=np.zeros(grid)
    rri=[0.0]*grid
    rrri=np.zeros(grid)
    em=np.zeros(grid)
    emom=np.zeros(grid)
    ne=grid-1
    wr=0
    
    for i in range(ne):
        if i==0:
            
            Pi[i]=P3
            rhoi[i]=rho3
            Ti[i]=T3
            hi[i]=h3
            si[i]=s3
            ui[i]=u3
            
            xii[i]=xi
            xif[i]=xf
            
            Ji[i]=J3
            Ni[i]=N3
            rcriti[i]=rcrit3
            Si[i]=Pi[i]/(p_sat(Ti[i],f,EoS))
            xgss[0]=1.001*P3
            xgss[1]=1.01*T3
            
        if i>0:
            xgss[0]=1.001*Pi[i]
            xgss[1]=1.01*Ti[i]
            
            
        sol=fsolve(Young_M,xgss,(Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,EoS,i),xtol=1.5e-015)
        diff=Young_M(sol,Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,EoS,i)
        print(diff)
        res=Young_res(sol,Pi[i],rhoi[i],hi[i],ui[i],rcriti[i],Ni,ri,xii[i],xif[i],f,EoS,i)
        em[i+1]=diff[0]
        emom[i+1]=diff[1]
        Pi[i+1]=res[0]
        rhoi[i+1]=res[1]
        Ti[i+1]=res[2]
        Si[i+1]=Pi[i+1]/(p_sat(Ti[i+1],f,EoS))
        hi[i+1]=res[3]
        si[i+1]=res[4]
        ui[i+1]=res[5]
        Yi[i+1]=res[6]
        delri[i+1]=res[7]
        rmi[i+1]=res[9]
        
        
        ri=res[8]
        rri[i+1]=ri
        nrate=Ncl_rate(Pi[i+1],Ti[i+1],f,EoS,k)
        Ji[i+1]=nrate[0]
        
        if wr==0:
            if Ji[i+1]<Ji[i]:
                print("shit")
                wr=1
                xc=xii[i]
        if wr==1:
            rnew[i+1]=rmi[i+1]
        rcriti[i+1]=nrate[1]
        
        xii[i+1]=xif[i]
        xif[i+1]=xii[i+1]+stp
        Ni[i+1]=Ji[i+1]*(xif[i+1]-xii[i+1])/(rhoi[i+1]*ui[i+1])/1000
        
        if i>1:
            rrri[i+1]=rri[i+1][1]
    """
    plt.figure(1)       
    plt.xlabel("Nozzle Axis")  
    plt.ylabel("Pressure (bar)")      
    plt.plot(xii,Pi,"--")
    plt.grid()
    """
    
    """
    plt.figure(2)       
    plt.xlabel("Nozzle Axis")  
    plt.ylabel("Pressure (bar)")      
    plt.plot(xii,Ji,"--")
    plt.grid()
    
    plt.figure(3)       
    plt.xlabel("Nozzle Axis")  
    plt.ylabel("Pressure (bar)")      
    plt.plot(xii,rmi,"--")
    plt.grid()
    """
    
    return ne,xii,Pi,Ti,ui,Yi,si,Si,Ji,rmi,rnew,em,emom,xc

##################################################################################################################

###################################################################################################################    




"""
#18c saturation
f="Water"
EoS="HEOS"
#P3=5832718.423752329
#T3=545.0486591928726
#u3=530.1670479959149
#xcond=0.9749999999999998


P3=5379618.3842589315
T3=535.7844060663691
u3=558.6371629074617
xcond=1.8749999999999993


Po=100.7*100000
A2=1
k=1

grid=1000
drpl=Gyarmathy_growth(P3,T3,u3,A2,grid,xcond,f,EoS,k)
"""


"""
exp=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Higpress_drpl/res18b/exp18b.csv")
plt.xlim(1,16)
plt.ylim(0.2,0.6)
plt.plot(exp["x"],exp["ppo"],"o")
plt.plot(drpl[1],drpl[2]/Po,"--")
"""




  
     
     
     
     
     
     
     

























    

"""
f="Water"
EoS="HEOS"

k=0.85

P3=9184.243731075647
T3=302.65355960937546
u3=561.3149166543993

J3=Ncl_rate(P3,T3,f,EoS,k)

rcrit3=J3[1]
J3=J3[0]



prop3=Cool_PT(P3,T3,f,EoS)
s3=prop3[4]
h3=prop3[3]
rho3=prop3[1]

xcond=50
stp=1/4

xi=xcond
xf=xi+stp

xgs=np.zeros(2)
xgs[0]=1.001*P3
xgs[1]=1.001*T3

P4=xgs[0]
T4=xgs[1]

"Gas properties"
propg=Cool_PT(P4,T4,f,EoS)
P=propg[0]
rhog=propg[1]
Tg=propg[2]
hg=propg[3]
sg=propg[4]
gamma=propg[7]
cpg=propg[9]
lamda=propg[10]
mug=propg[11]
Prg=cpg*mug/lamda


propsat_p=Cool_Psat(P4,f,EoS)
Tsat=propsat_p[1]

rhol=1/propsat_p[6]
hl=propsat_p[7]
sl=propsat_p[8]
    
    
L=hg-hl
R=8.31462/CP.PropsSI('molar_mass',f)

"Droplet growth"
alpha=0
theta=R*Tsat/L*(alpha-2-1*(gamma+1)/(2*gamma)*(cpg*Tsat/L))
delT=Tsat-Tg

A=lamda*delT
B=3.78*(1-theta)/Prg
C=rhol*L
D=1.5*mug*math.sqrt(R*Tg)/(2*P)
E=rcrit3
F=(xf-xi)/u3


coeff=np.zeros(4)
coeff[0]=C
coeff[1]=C*B*D-C*E
coeff[2]=-(C*B*D*E+F*A)
coeff[3]=F*A*rcrit3

roots=np.roots(coeff)
rm=max(roots)

delr=F*A/C*(1-rcrit3/rm)/(rm+B*D)
r1=rcrit3+delr
"""








