# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 12:39:44 2019

@author: USER
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from Twophase import fnd_2phase_exp, Propmix, sol2Phase
from Coollib import Cool_PT, Cool_Ps, Cool_hs, Cool_Psat


def Rayleigh(x,P3,rho3,u3,P4,rhog):
    rho4=rhog/x
    y=P4+(rho3*u3)**2.0/rho4-P3-rho3*u3**2.0
    return y

def Rayleigh_rho(P3,rho3,u3,P4,rhog):
    guessx=0.001
    ne=100000
    stpx=1/250
    
    err=np.zeros(ne)
    
    for i in range (ne):
        err[i]=Rayleigh(guessx,P3,rho3,u3,P4,rhog)
        if i>0:
            if err[i]*err[i-1]<0:
                guessx=guessx-0.5*stpx
                X=fsolve(Rayleigh,guessx,(P3,rho3,u3,P4,rhog))
                errorx=Rayleigh(X,P3,rho3,u3,P4,rhog)
                print (X)
                if abs (errorx>1.8e-08):
                    print ("Quality not found for",P4,errorx)
                           
                break
        
        guessx=guessx+stpx
        
    rho4=rhog/X
    return rho4



def Rayleighl(P3,T3,u3,f,EoS):
    prop3=Cool_PT(P3,T3,f,EoS)
    rho3=prop3[1]
    
    
    Pf=2.75*P3
    ne=1000
    stp=(Pf-P3)/ne
    P4=P3+stp
    
    Pi=np.zeros(ne)
    Vi=np.zeros(ne)
    
    for i in range (ne):
        
        prop4=Cool_Psat(P4,f,EoS)
        rhog=1/prop4[2]
        rho4=Rayleigh_rho(P3,rho3,u3,P4,rhog)
        Pi[i]=P4
        Vi[i]=1/rho4
        
        print(i)
        P4=P4+stp
        
        
        

    return Vi,Pi

def diff_Hugoniot(X4,P3,rho3,h3,P4,f,EoS,q):
    prop4=Propmix(P4,X4,f,EoS)
    h4=prop4[3]
    rho4=prop4[1]
    y=-q+h4-h3-0.5*(1/rho3+1/rho4)*(P4-P3)
    return y

    

def fnd_Hugoniot(P4,P3,T3,q,f,EoS):    


    prop3=Cool_PT(P3,T3,f,EoS)
    rho3=prop3[1]
    h3=prop3[3]
    
       
    xi=0.4
    xf=1.5
    ne=200
    stpx=(xf-xi)/ne
    
    err=np.zeros(ne)
    for i in range (ne):
        err[i]=diff_Hugoniot(xi,P3,rho3,h3,P4,f,EoS,q)
        if i>0:
            if (err[i]*err[i-1])<0:
                xi=xi-0.5*stpx
                X4=fsolve(diff_Hugoniot,xi,(P3,rho3,h3,P4,f,EoS,q))
                error=diff_Hugoniot(X4,P3,rho3,h3,P4,f,EoS,q)
                if error>1.8e-07:
                    print ("hugoniot no found for",P4,error)
                break
        xi=xi+stpx
    
    prop4=Propmix(P4,X4,f,EoS)
    rho4=prop4[1]
    return rho4


def Hugoniotl(P3,T3,q,f,EoS):

    Pf=3*P3
    ne=100
    stp=(Pf-P3)/ne
    P4=P3+stp
    
    Pi=np.zeros(ne)
    Vi=np.zeros(ne)
    
    for i in range (ne):
        rhoi=fnd_Hugoniot(P4,P3,T3,q,f,EoS)
        Vi[i]=1/rhoi
        Pi[i]=P4/100000
        
        P4=P4+stp
    
    
    

    return Vi,Pi

def diffRH_EQ(P4,P3,T3,u3,q,f,EoS):
    prop3=Cool_PT(P3,T3,f,EoS)
    rho3=prop3[1]
    h3=prop3[3]
    
    prop4=Cool_Psat(P4,f,EoS)
    rhog=1/prop4[2]
    
    X4=(P3-P4+rho3*u3**2.0)*rhog/(rho3*u3)**2.0
    prop4=Propmix(P4,X4,f,EoS)
    h4=prop4[3]
    y=-q+h4-h3-0.5*(1/rho3+X4/rhog)*(P4-P3)
    return y

def fndX_RHEQ(P4,P3,T3,u3,f,EoS):
  
    prop3=Cool_PT(P3,T3,f,EoS)
    rho3=prop3[1]
    prop4=Cool_Psat(P4,f,EoS)
    rhog=1/prop4[2]
    
    X4=(P3-P4+rho3*u3**2.0)*rhog/(rho3*u3)**2.0
    return X4

    

def fndPs_RHeq(P3,T3,u3,q,f,EoS):
    
    pi=1*P3
    pf=3*P3
    ne=5000
    stp=(pf-pi)/ne
    err=np.zeros(ne)
    
    
    for i in range (ne):
        err[i]=diffRH_EQ(pi,P3,T3,u3,q,f,EoS)
        if i==(ne-1):
            Pwd=0
        if i >0:
            if (err[i]*err[i-1])<0:
                gssp=pi-0.5*stp
                Pwd=fsolve(diffRH_EQ,gssp,(P3,T3,u3,q,f,EoS))
                error=diffRH_EQ(Pwd,P3,T3,u3,q,f,EoS)
                if error>1.8e-07:
                    print("error Pwd not found", error)
                
                break
        
        
        pi=pi+stp
        
    
    for i in range (ne):
        err[i]=diffRH_EQ(pi,P3,T3,u3,q,f,EoS)
        if i==(ne-1):
            Psd=0
        if i >0:
            if (err[i]*err[i-1])<0:
                gssp=pi-0.5*stp
                Psd=fsolve(diffRH_EQ,gssp,(P3,T3,u3,q,f,EoS))
                error=diffRH_EQ(Psd,P3,T3,u3,q,f,EoS)
                if error>1.8e-07:
                    print("error Psd not found", error)
                
                break
        pi=pi+stp
        
    return Pwd, Psd
   

def CJi(q,stp,P3,T3,u3,f,EoS):
    qCJ=0
    ne=10000000 
    for i in range (ne):
        Ps=fndPs_RHeq(P3,T3,u3,q,f,EoS)
        Pwd=Ps[0]
        Psd=Ps[1]
        diff=Pwd-Psd
        print (diff)
        if Ps[0]==0:
            qnew=q-stp
            stpn=stp/10
            break
        if abs(diff)<10.92:
            if abs(diff)>0:
                
                qnew=q
                stpn=stp
                qCJ=q
                break
        
        q=q+stp
        
    return qnew,stpn,qCJ,diff



def fnd_CJ(P3,T3,u3,f,EoS):
    q=0
    stp=10000
    
    for i in range(10000000000):
        sol=CJi(q,stp,P3,T3,u3,f,EoS)
        if sol[2]!=0:
            qCJ=sol[2]
            break
        
        q=sol[0]
        stp=sol[1]
        print (stp,q,sol[3])
        
    Ps=fndPs_RHeq(P3,T3,u3,qCJ,f,EoS)
    Pwd_CJ=Ps[0]
    Psd_CJ=Ps[1]
    Xwd_CJ=fndX_RHEQ(Pwd_CJ,P3,T3,u3,f,EoS)
    Xsd_CJ=fndX_RHEQ(Psd_CJ,P3,T3,u3,f,EoS)
    wd_CJ=Propmix(Pwd_CJ,Xwd_CJ,f,EoS)
    sd_CJ=Propmix(Psd_CJ,Xsd_CJ,f,EoS)
    Vwd=1/wd_CJ[1]
    Vsd=1/sd_CJ[1]
    swd=wd_CJ[4]
    ssd=sd_CJ[4]
    diffV=Vsd-Vwd
    return Pwd_CJ,Xwd_CJ,Vwd,swd,qCJ,Vsd,Psd_CJ,diffV

def diff_RayleighS(guess_X,P4,P3,rho3,J3_2,f,EoS):
    
    prop4=Propmix(P4,guess_X,f,EoS)
    rho4=prop4[1]
    
    y=-J3_2+(P4-P3)/(1/rho3-1/rho4)
    return y


def fnd_RayleighS(P4,P3,rho3,J3_2,f,EoS):
    

    guess_Xi=0.01
    guess_Xf=2.5
    ne=1000
    stp=(guess_Xf-guess_Xi)/ne
    errR=np.zeros(ne)
    
    
    for i in range(ne):
        errR[i]=diff_RayleighS(guess_Xi,P4,P3,rho3,J3_2,f,EoS)
        if i>0:
            if (errR[i]*errR[i-1])<0:
                guess_Xi=guess_Xi-0.5*stp
                X4=fsolve(diff_RayleighS,guess_Xi,(P4,P3,rho3,J3_2,f,EoS))
                break
            
        guess_Xi=guess_Xi+stp
    
    Prop4=Propmix(P4,X4,f,EoS)
    rho4=Prop4[1]
  
        
    
    return X4,rho4



    
def RayleighlS(qR,P3,T3,u3,f,EoS):
    
    Prop3=Cool_PT(P3,T3,f,EoS)
    rho3=Prop3[1]
    J3=qR*rho3*u3
    J3_2=J3*J3
    
    
    Pi=P3+0.01
    Pf=2.5*P3
    ne=50
    stp=(Pf-Pi)/ne
    
    P=np.zeros(ne)
    V=np.zeros(ne)
    X=np.zeros(ne)
    
    Pi=1.01*P3
    for i in range (ne):
        R4=fnd_RayleighS(Pi,P3,rho3,J3_2,f,EoS)
        P[i]=Pi/100000
        X[i]=R4[0]
        V[i]=1/R4[1]
        
        print (i)
        Pi=Pi+stp
    
    return V,P


def diff_RHS(guessP4,P3,T3,u3,qR,f,EoS):
    Prop3=Cool_PT(P3,T3,f,EoS)
    rho3=Prop3[1]
    h3=Prop3[3]
    J3=qR*rho3*u3
    J3_2=J3*J3
    
    
    
    
    P4=guessP4
    prop4=Cool_Psat(P4,f,EoS)
    rhog=1/prop4[2]
    X4=((P3-P4)/(qR*J3_2)+1/rho3)*rhog
    prop4=Propmix(P4,X4,f,EoS)
    h4=prop4[3]
    y=h4-h3-0.5*(1/rho3+X4/rhog)*(P4-P3)
    return y        


def fndX_RHS(guessP4,P3,T3,u3,qR,f,EoS):
    Prop3=Cool_PT(P3,T3,f,EoS)
    rho3=Prop3[1]
    J3=qR*rho3*u3
    J3_2=J3*J3
    
        
    P4=guessP4
    prop4=Cool_Psat(P4,f,EoS)
    rhog=1/prop4[2]
    X4=((P3-P4)/(qR*J3_2)+1/rho3)*rhog
    
    return X4






def fnd_RHS(P3,T3,u3,qR,f,EoS):
    

    guessP4i=1*P3
    guessP5f=3*P3
    ne=7000
    stp=(guessP5f-guessP4i)/ne
    y1=np.zeros(ne)
    y2=np.zeros(ne)
    
    for i in range (ne):
       y1[i]=diff_RHS(guessP4i,P3,T3,u3,qR,f,EoS) 
       if i>0:
           if (y1[i]*y1[i-1])<0:
               guessP4i=guessP4i-0.5*stp
               Pwd=fsolve(diff_RHS,guessP4i,(P3,T3,u3,qR,f,EoS))
               err=diff_RHS(Pwd,P3,T3,u3,qR,f,EoS)
               if err>1.7e-8:
                   print ("Pwd not found", err)
               
               Xwd=fndX_RHS(Pwd,P3,T3,u3,qR,f,EoS)
               propwd=Propmix(Pwd,Xwd,f,EoS)
               Vwd=1/propwd[1]
               guessP4i=guessP4i+stp
               break
           if i==(ne-1):
               Pwd=0
               Vwd=0
               Xwd=0
               
       guessP4i=guessP4i+stp
    
    for i in range (ne):
       y2[i]=diff_RHS(guessP4i,P3,T3,u3,qR,f,EoS) 
       if i>0:
           if (y2[i]*y2[i-1])<0:
               guessP4i=guessP4i-0.5*stp
               Psd=fsolve(diff_RHS,guessP4i,(P3,T3,u3,qR,f,EoS))
               Xsd=fndX_RHS(Psd,P3,T3,u3,qR,f,EoS)
               propsd=Propmix(Psd,Xsd,f,EoS)
               Vsd=1/propsd[1]
               break
           if i==(ne-1):
               Psd=0
               Xsd=0
               Vsd=0
               
           
        
       guessP4i=guessP4i+stp
       
     
       
    return Vwd,Pwd,Xwd,Vsd,Psd,Xsd
       

def CJS(gssPR,stpi,P3,T3,u3,f,EoS):
    j=0
    lol1=fnd_RHS(P3,T3,u3,gssPR,f,EoS)
    for i in range (1000):
        lol1=fnd_RHS(P3,T3,u3,gssPR,f,EoS)
        Vwd=lol1[0]
        Vsd=lol1[3]
        Pwd=lol1[1]
        Psd=lol1[4]
        d=abs(Pwd-Psd)
        print (d,gssPR)
        if d==0:
            gssPR=gssPR-stpi
            stpi=stpi/10
            VCJ=0
            PCJ=0
            qR=0
            break
        
        if abs (d)<20:
            if d>0:
                VCJ=Vwd
                PCJ=lol1[1]
                qR=gssPR
                j=1
                break
                
                
        
            
        gssPR=gssPR+stpi
    return gssPR, stpi, VCJ, PCJ,d,qR,j


def fnd_CJS(P3,T3, u3, f, EoS):
    
    gssPR=1
    stpi=-0.1
    
    for i in range (20):
        cj=CJS(gssPR,stpi,P3,T3,u3,f,EoS)
        gssPR=cj[0]
        if abs( cj[6])==1:
            VCJ=cj[2]
            PCJ=cj[3]
            JPR=cj[5]
            break
        gssPR=cj[0]
        stpi=cj[1]
        print(i)
        
    return VCJ, PCJ, JPR






"""
  
f="Water"
EoS="PR"

P3=29537.218652288513
T3=305.3151097579453
u3=556.2981851531378

Prop3=Cool_PT(P3,T3,f,EoS)
V3=1/Prop3[1]

#cj=fnd_CJS(P3,T3,u3,f,EoS)
Vcj=3.717845795071342
Pcj=0.4301208497335008
Jcj=0.9836204900000003

det=fnd_RHS(P3,T3,u3,1,f,EoS)
Vwd=det[0]
Pwd=det[1]/100000

Vsd=det[3]
Psd=det[4]/100000



cj_c=fnd_RHS(P3,T3,u3,Jcj,f,EoS)
Vcj_c=float(cj_c[0])
Pcj_c=float(cj_c[1])/100000
print(Vcj_c-Vcj,Pcj_c-Pcj)
rl=RayleighlS(1,P3,T3,u3,f,EoS)
rlcj=RayleighlS(Jcj-0.01,P3,T3,u3,f,EoS)
hlcj=Hugoniotl(P3,T3,0,f,EoS)
P3=P3/100000


plt.figure(figsize=(7,5))

family="serif"
axisize=14
lsize=16



font = {'family': family,
        'weight': 'normal',
        'size': 12,
        }

plt.ylabel("Pressure (bar)",fontsize=lsize,family=family)
plt.xlabel("Specific volume (m$^{3}$/kg)",fontsize=lsize,family=family)
plt.ylim(0.255, 0.6)
plt.xlim(3, 5)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)



plt.plot(hlcj[0],hlcj[1],"-",color="darkcyan",lw=2,label="Hugoniot")

plt.plot(rl[0],rl[1],"-",color="darkred",lw=1.5,label="Rayleigh")
plt.plot(rlcj[0],rlcj[1],"-",color="darkred",lw=1.5)



plt.plot(V3,P3,"v",ms=8,color="indigo",label="State 3")
plt.plot(Vwd,Pwd,"o",ms=7,color="chocolate",label="Weak detonation")
plt.plot(Vsd,Psd,"o",ms=7,color="darkgreen",label="Strong detonation")
plt.plot(Vcj,Pcj,"^",ms=8,color="black",label="CJ-point")
plt.legend(frameon=False,prop=font)
plt.text(3.05,0.445,"J=0.983J$_{3}$",fontsize=12,family=family)
plt.text(3.55,0.5,"J=J$_{3}$",fontsize=12,family=family)
plt.arrow(3.6,0.495,0,-0.03,head_width=0.03, head_length=0.01,fc='k', ec='k')
plt.arrow(3.225,0.46,0,+0.02,head_width=0.032, head_length=0.012,fc='k', ec='k')
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.savefig('RH1.pdf')             
"""

###########################################################################################
###########################################################################################
###IF95-B107.4-RH
"""
f="Water"
EoS="HEOS"

P3=29936.74807487126
T3=309.2909164923855
u3=547.0851645209306


P3=29537.218652288513
T3=305.3151097579453
u3=556.2981851531378

Prop3=Cool_PT(P3,T3,f,EoS)
V3=1/Prop3[1]

#cj=fnd_CJS(P3,T3,u3,f,EoS)
Vcj=3.650863323886017753e+00
Pcj=4.277717599011196580e+04/100000
Jcj=0.9743370000000002

det=fnd_RHS(P3,T3,u3,1,f,EoS)
Vwd=det[0]
Pwd=det[1]/100000

Vsd=det[3]
Psd=det[4]/100000



cj_c=fnd_RHS(P3,T3,u3,Jcj,f,EoS)
Vcj_c=float(cj_c[0])
Pcj_c=float(cj_c[1])/100000
print(Vcj_c-Vcj,Pcj_c-Pcj)
rl=RayleighlS(1,P3,T3,u3,f,EoS)
rlcj=RayleighlS(Jcj-0.015,P3,T3,u3,f,EoS)
hlcj=Hugoniotl(P3,T3,0,f,EoS)
P3=P3/100000


plt.figure(figsize=(7,5))

family="serif"
axisize=14
lsize=16



font = {'family': family,
        'weight': 'normal',
        'size': 12,
        }

plt.ylabel("Pressure (bar)",fontsize=lsize,family=family)
plt.xlabel("Specific volume (m$^{3}$/kg)",fontsize=lsize,family=family)
plt.ylim(0.255, 0.6)
plt.xlim(3, 5)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)



plt.plot(hlcj[0],hlcj[1],"-",color="darkcyan",lw=2,label="Hugoniot")

plt.plot(rl[0],rl[1],"-",color="darkred",lw=1.5,label="Rayleigh")
plt.plot(rlcj[0],rlcj[1],"-",color="darkred",lw=1.5)



plt.plot(V3,P3,"v",ms=8,color="cornflowerblue",label="State 3")
plt.plot(Vwd,Pwd,"o",ms=7,color="chocolate",label="Weak detonation")
plt.plot(Vsd,Psd,"o",ms=7,color="darkgreen",label="Strong detonation")
plt.plot(Vcj,Pcj,"^",ms=8,color="black",label="CJ-point")
plt.legend(frameon=False,prop=font)
plt.text(3.05,0.43,"J=0.9743J$_{3}$",fontsize=12,family=family)
plt.text(3.55,0.5,"J=J$_{3}$",fontsize=12,family=family)
plt.arrow(3.6,0.490,0,-0.03,head_width=0.03, head_length=0.01,fc='k', ec='k')
plt.arrow(3.225,0.45,0,+0.02,head_width=0.032, head_length=0.012,fc='k', ec='k')
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.savefig('RHIF951074.pdf')
"""
###################################################################################
###################################################################################
###IF95-112.2

"""
f="Water"
EoS="HEOS"

P3=27380.717149493506
T3=307.1414656148738
u3=572.7001381579616

Prop3=Cool_PT(P3,T3,f,EoS)
V3=1/Prop3[1]

#cj=fnd_CJS(P3,T3,u3,f,EoS)
Vcj=4.03063986
Pcj=38809.65961684/100000
Jcj=0.9494402000000004

det=fnd_RHS(P3,T3,u3,1,f,EoS)
Vwd=det[0]
Pwd=det[1]/100000

Vsd=det[3]
Psd=det[4]/100000



cj_c=fnd_RHS(P3,T3,u3,Jcj,f,EoS)
Vcj_c=float(cj_c[0])
Pcj_c=float(cj_c[1])/100000
print(Vcj_c-Vcj,Pcj_c-Pcj)
rl=RayleighlS(1,P3,T3,u3,f,EoS)
rlcj=RayleighlS(Jcj-0.025,P3,T3,u3,f,EoS)
hlcj=Hugoniotl(P3,T3,0,f,EoS)
P3=P3/100000


plt.figure(figsize=(7,5))

family="serif"
axisize=14
lsize=16



font = {'family': family,
        'weight': 'normal',
        'size': 12,
        }

plt.ylabel("Pressure (bar)",fontsize=lsize,family=family)
plt.xlabel("Specific volume (m$^{3}$/kg)",fontsize=lsize,family=family)
plt.ylim(0.255, 0.6)
plt.xlim(3, 5.25)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)



plt.plot(hlcj[0],hlcj[1],"-",color="darkcyan",lw=2,label="Hugoniot")

plt.plot(rl[0],rl[1],"-",color="darkred",lw=1.5,label="Rayleigh")
plt.plot(rlcj[0],rlcj[1],"-",color="darkred",lw=1.5)



plt.plot(V3,P3,"v",ms=8,color="limegreen",label="State 3")
plt.plot(Vwd,Pwd,"o",ms=7,color="chocolate",label="Weak detonation")
plt.plot(Vsd,Psd,"o",ms=7,color="darkgreen",label="Strong detonation")
plt.plot(Vcj,Pcj,"^",ms=8,color="black",label="CJ-point")
plt.legend(frameon=False,prop=font)
plt.text(3.05,0.42,"J=0.9494J$_{3}$",fontsize=12,family=family)
plt.text(3.55,0.51,"J=J$_{3}$",fontsize=12,family=family)
plt.arrow(3.6,0.505,0,-0.03,head_width=0.03, head_length=0.01,fc='k', ec='k')
plt.arrow(3.225,0.44,0,+0.02,head_width=0.032, head_length=0.012,fc='k', ec='k')
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.savefig('RHIF951122.pdf')
"""
###################################################################################
###################################################################################
###IF95-115.2
"""
f="Water"
EoS="HEOS"

P3=25888.55185297716
T3=305.8005465467396
u3=588.2523896809308

Prop3=Cool_PT(P3,T3,f,EoS)
V3=1/Prop3[1]
#cj=fnd_CJS(P3,T3,u3,f,EoS)
Vcj=4.24167725
Pcj=0.36738542
Jcj=0.93176358

det=fnd_RHS(P3,T3,u3,1,f,EoS)
Vwd=det[0]
Pwd=det[1]/100000

Vsd=det[3]
Psd=det[4]/100000



cj_c=fnd_RHS(P3,T3,u3,Jcj,f,EoS)
Vcj_c=float(cj_c[0])
Pcj_c=float(cj_c[1])/100000
print(Vcj_c-Vcj,Pcj_c-Pcj)
rl=RayleighlS(1,P3,T3,u3,f,EoS)
rlcj=RayleighlS(Jcj-0.032,P3,T3,u3,f,EoS)
hlcj=Hugoniotl(P3,T3,0,f,EoS)
P3=P3/100000


plt.figure(figsize=(7,5))

family="serif"
axisize=14
lsize=16



font = {'family': family,
        'weight': 'normal',
        'size': 12,
        }

plt.ylabel("Pressure (bar)",fontsize=lsize,family=family)
plt.xlabel("Specific volume (m$^{3}$/kg)",fontsize=lsize,family=family)
plt.ylim(0.225, 0.6)
plt.xlim(2.75, 5.5)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)



plt.plot(hlcj[0],hlcj[1],"-",color="darkcyan",lw=2,label="Hugoniot")

plt.plot(rl[0],rl[1],"-",color="darkred",lw=1.5,label="Rayleigh")
plt.plot(rlcj[0],rlcj[1],"-",color="darkred",lw=1.5)



plt.plot(V3,P3,"v",ms=8,color="crimson",label="State 3")
plt.plot(Vwd,Pwd,"o",ms=7,color="chocolate",label="Weak detonation")
plt.plot(Vsd,Psd,"o",ms=7,color="darkgreen",label="Strong detonation")
plt.plot(Vcj,Pcj,"^",ms=8,color="black",label="CJ-point")
plt.legend(frameon=False,prop=font)
plt.text(3.05,0.41,"J=0.9317J$_{3}$",fontsize=12,family=family)
plt.text(3.55,0.52,"J=J$_{3}$",fontsize=12,family=family)
plt.arrow(3.6,0.51,0,-0.03,head_width=0.03, head_length=0.01,fc='k', ec='k')
plt.arrow(3.225,0.43,0,+0.02,head_width=0.032, head_length=0.012,fc='k', ec='k')
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.savefig('RHIF951152.pdf')
"""
###################################################################################
###################################################################################
###IF95-119.7
"""
f="Water"
EoS="HEOS"

P3=23778.859569940047
T3=303.77303607841304
u3=611.1138930246722

Prop3=Cool_PT(P3,T3,f,EoS)
V3=1/Prop3[1]

#cj=fnd_CJS(P3,T3,u3,f,EoS)
Vcj=4.58210015
Pcj=33811.79247113/100000
Jcj=0.9071155900000004


det=fnd_RHS(P3,T3,u3,1,f,EoS)
Vwd=det[0]
Pwd=det[1]/100000

Vsd=det[3]
Psd=det[4]/100000



cj_c=fnd_RHS(P3,T3,u3,Jcj,f,EoS)
Vcj_c=float(cj_c[0])
Pcj_c=float(cj_c[1])/100000
print(Vcj_c-Vcj,Pcj_c-Pcj)
rl=RayleighlS(1,P3,T3,u3,f,EoS)
rlcj=RayleighlS(Jcj-0.045,P3,T3,u3,f,EoS)
hlcj=Hugoniotl(P3,T3,0,f,EoS)
P3=P3/100000


plt.figure(figsize=(7,5))

family="serif"
axisize=14
lsize=16



font = {'family': family,
        'weight': 'normal',
        'size': 12,
        }

plt.ylabel("Pressure (bar)",fontsize=lsize,family=family)
plt.xlabel("Specific volume (m$^{3}$/kg)",fontsize=lsize,family=family)
plt.ylim(0.225, 0.6)
plt.xlim(2.75, 6)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)



plt.plot(hlcj[0],hlcj[1],"-",color="darkcyan",lw=2,label="Hugoniot")

plt.plot(rl[0],rl[1],"-",color="darkred",lw=1.5,label="Rayleigh")
plt.plot(rlcj[0],rlcj[1],"-",color="darkred",lw=1.5)



plt.plot(V3,P3,"v",ms=8,color="dodgerblue",label="State 3")
plt.plot(Vwd,Pwd,"o",ms=7,color="chocolate",label="Weak detonation")
plt.plot(Vsd,Psd,"o",ms=7,color="darkgreen",label="Strong detonation")
plt.plot(Vcj,Pcj,"^",ms=8,color="black",label="CJ-point")
plt.legend(frameon=False,prop=font)
plt.text(3.05,0.40,"J=0.9071J$_{3}$",fontsize=12,family=family)
plt.text(3.55,0.535,"J=J$_{3}$",fontsize=12,family=family)
plt.arrow(3.6,0.525,0,-0.03,head_width=0.03, head_length=0.01,fc='k', ec='k')
plt.arrow(3.225,0.42,0,+0.02,head_width=0.032, head_length=0.012,fc='k', ec='k')
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.savefig('RHIF951197.pdf')
"""
######################################################################################
#####################################################################################
###PR-119.7
"""
f="Water"
EoS="PR"

P3=27604.82760978155
T3=303.7792529538604
u3=576.3613466359263

Prop3=Cool_PT(P3,T3,f,EoS)
V3=1/Prop3[1]
#cj=fnd_CJS(P3,T3,u3,f,EoS)
Vcj=3.955988
Pcj=40253.07884453/100000
Jcj=0.9596373000000006


det=fnd_RHS(P3,T3,u3,1,f,EoS)
Vwd=det[0]
Pwd=det[1]/100000

Vsd=det[3]
Psd=det[4]/100000



cj_c=fnd_RHS(P3,T3,u3,Jcj,f,EoS)
Vcj_c=float(cj_c[0])
Pcj_c=float(cj_c[1])/100000
print(Vcj_c-Vcj,Pcj_c-Pcj)
rl=RayleighlS(1,P3,T3,u3,f,EoS)
rlcj=RayleighlS(Jcj-0.02,P3,T3,u3,f,EoS)
hlcj=Hugoniotl(P3,T3,0,f,EoS)
P3=P3/100000


plt.figure(figsize=(7,5))

family="serif"
axisize=14
lsize=16



font = {'family': family,
        'weight': 'normal',
        'size': 12,
        }

plt.ylabel("Pressure (bar)",fontsize=lsize,family=family)
plt.xlabel("Specific volume (m$^{3}$/kg)",fontsize=lsize,family=family)
plt.ylim(0.225, 0.6)
plt.xlim(2.75, 5.25)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)



plt.plot(hlcj[0],hlcj[1],"-",color="darkcyan",lw=2,label="Hugoniot")

plt.plot(rl[0],rl[1],"-",color="darkred",lw=1.5,label="Rayleigh")
plt.plot(rlcj[0],rlcj[1],"-",color="darkred",lw=1.5)



plt.plot(V3,P3,"v",ms=8,color="dodgerblue",label="State 3")
plt.plot(Vwd,Pwd,"^",ms=8,color="maroon",label="Weak detonation")
plt.plot(Vsd,Psd,"^",ms=8,color="midnightblue",label="Strong detonation")
plt.plot(Vcj,Pcj,"o",ms=6,color="black",label="CJ-point")
plt.legend(frameon=False,prop=font)
plt.text(3.05,0.43,"J=0.9596J$_{3}$",fontsize=12,family=family)
plt.text(3.55,0.52,"J=J$_{3}$",fontsize=12,family=family)
plt.arrow(3.6,0.51,0,-0.03,head_width=0.03, head_length=0.01,fc='k', ec='k')
plt.arrow(3.225,0.452,0,+0.02,head_width=0.032, head_length=0.012,fc='k', ec='k')
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.savefig('RHPR1197.pdf')
"""
########################################################################################
##########################################################################################
###PR-115.2

"""

f="Water"
EoS="PR"

P3=29569.543042423236
T3=305.3979487151169
u3=556.022208961323

Prop3=Cool_PT(P3,T3,f,EoS)
V3=1/Prop3[1]
#cj=fnd_CJS(P3,T3,u3,f,EoS)
Vcj=3.71522642
Pcj=43046.08200542/100000
Jcj=0.9838574000000004

det=fnd_RHS(P3,T3,u3,1,f,EoS)
Vwd=det[0]
Pwd=det[1]/100000

Vsd=det[3]
Psd=det[4]/100000



cj_c=fnd_RHS(P3,T3,u3,Jcj,f,EoS)
Vcj_c=float(cj_c[0])
Pcj_c=float(cj_c[1])/100000
print(Vcj_c-Vcj,Pcj_c-Pcj)
rl=RayleighlS(1,P3,T3,u3,f,EoS)
rlcj=RayleighlS(Jcj-0.01,P3,T3,u3,f,EoS)
hlcj=Hugoniotl(P3,T3,0,f,EoS)
P3=P3/100000


plt.figure(figsize=(7,5))

family="serif"
axisize=14
lsize=16



font = {'family': family,
        'weight': 'normal',
        'size': 12,
        }

plt.ylabel("Pressure (bar)",fontsize=lsize,family=family)
plt.xlabel("Specific volume (m$^{3}$/kg)",fontsize=lsize,family=family)
plt.ylim(0.275, 0.55)
plt.xlim(3, 4.8)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)



plt.plot(hlcj[0],hlcj[1],"-",color="darkcyan",lw=2,label="Hugoniot")

plt.plot(rl[0],rl[1],"-",color="darkred",lw=1.5,label="Rayleigh")
plt.plot(rlcj[0],rlcj[1],"-",color="darkred",lw=1.5)



plt.plot(V3,P3,"v",ms=8,color="crimson",label="State 3")
plt.plot(Vwd,Pwd,"^",ms=8,color="maroon",label="Weak detonation")
plt.plot(Vsd,Psd,"^",ms=8,color="midnightblue",label="Strong detonation")
plt.plot(Vcj,Pcj,"o",ms=6,color="black",label="CJ-point")
plt.legend(frameon=False,prop=font)
plt.text(3.05,0.44,"J=0.9838J$_{3}$",fontsize=12,family=family)
plt.text(3.55,0.50,"J=J$_{3}$",fontsize=12,family=family)
plt.arrow(3.6,0.498,0,-0.03,head_width=0.03, head_length=0.01,fc='k', ec='k')
plt.arrow(3.225,0.46,0,+0.02,head_width=0.032, head_length=0.012,fc='k', ec='k')
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.savefig('RHPR1152.pdf')

"""

####################################################################################
###PR-112.2
"""
f="Water"
EoS="PR"

P3=30974.59408384412
T3=306.5057643942014
u3=541.9549302290164

rl=RayleighlS(1,P3,T3,u3,f,EoS)
hlcj=Hugoniotl(P3,T3,0,f,EoS)

plt.figure(figsize=(7,5))

family="serif"
axisize=14
lsize=16



font = {'family': family,
        'weight': 'normal',
        'size': 12,
        }

plt.ylabel("Pressure (bar)",fontsize=lsize,family=family)
plt.xlabel("Specific volume (m$^{3}$/kg)",fontsize=lsize,family=family)
plt.ylim(0.44, 0.46)
plt.xlim(3.5, 3.65)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)
plt.plot(hlcj[0],hlcj[1],"-",color="darkcyan",lw=2,label="Hugoniot")
plt.plot(rl[0],rl[1],"-",color="darkred",lw=1.5,label="Rayleigh")
plt.legend(frameon=False,prop=font)
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.savefig('RHPR1122.pdf')
"""
########################################################################33
########PR107.4
"""
f="Water"
EoS="PR"

P3=33394.896987440174
T3=308.3246933849502
u3=518.5134033068946

rl=RayleighlS(1,P3,T3,u3,f,EoS)
hlcj=Hugoniotl(P3,T3,0,f,EoS)

plt.figure(figsize=(7,5))

family="serif"
axisize=14
lsize=16



font = {'family': family,
        'weight': 'normal',
        'size': 12,
        }

plt.ylabel("Pressure (bar)",fontsize=lsize,family=family)
plt.xlabel("Specific volume (m$^{3}$/kg)",fontsize=lsize,family=family)
plt.ylim(0.35, 0.55)
plt.xlim(3, 4.5)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)
plt.plot(hlcj[0],hlcj[1],"-",color="darkcyan",lw=2,label="Hugoniot")
plt.plot(rl[0],rl[1],"-",color="darkred",lw=1.5,label="Rayleigh")
plt.legend(frameon=False,prop=font)
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.savefig('RHPR1074.pdf')
"""








