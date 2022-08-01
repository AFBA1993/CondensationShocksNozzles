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
    Pf=2.8*P3
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

def RayleighlS1(qR,P3,T3,u3,f,EoS):
    
    Prop3=Cool_PT(P3,T3,f,EoS)
    rho3=Prop3[1]
    J3=qR*rho3*u3
    J3_2=J3*J3
    
    
    Pi=P3+0.01
    Pf=2*P3
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

#####################################################
############# Example ###############################
#####################################################

"""
f="Water"
EoS="HEOS"

P3=21114.352644768645
T3=303.3289479725017
u3=539.5920291743387
Prop3=Cool_PT(P3,T3,f,EoS)
V3=1/Prop3[1]
#cj=fnd_CJS(P3,T3,u3,f,EoS)
Vcj=5.20078796
Pcj=29649.64279269/100000
Jcj=0.9798407000000001

det=fnd_RHS(P3,T3,u3,1,f,EoS)
Vwd=det[0]
Pwd=det[1]/100000

Vsd=det[3]
Psd=det[4]/100000


rl=RayleighlS(1,P3,T3,u3,f,EoS)
rlcj=RayleighlS(Jcj-0.01,P3,T3,u3,f,EoS)
hlcj=Hugoniotl(P3,T3,0,f,EoS)
P3=P3/100000


plt.figure(figsize=(7.5,5))

family="serif"
axisize=14
lsize=16



font = {'family': family,
        'weight': 'normal',
        'size': 12,
        }

plt.ylabel("Pressure (bar)",fontsize=lsize,family=family)
plt.xlabel("Specific volume (m$^{3}$/kg)",fontsize=lsize,family=family)
plt.ylim(0.2, 0.4)
plt.xlim(4, 7)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)



plt.plot(hlcj[0],hlcj[1],"-",color="darkcyan",lw=2,label="Hugoniot")

plt.plot(rl[0],rl[1],"-",color="darkred",lw=1.5,label="Rayleigh")
plt.plot(rlcj[0],rlcj[1],"-",color="darkred",lw=1.5)



plt.plot(V3,P3,"v",ms=8,color="limegreen",label="Upstream state")
plt.plot(Vwd,Pwd,"o",ms=7,color="blue",label="Weak detonation")
plt.plot(Vsd,Psd,"o",ms=7,color="darkgreen",label="Strong detonation")
plt.plot(Vcj,Pcj,"^",ms=8,color="black",label="CJ-point")
plt.legend(frameon=False,prop=font)

plt.text(4.6,0.348,"\u03C1u=(\u03C1u)$_{US}$",fontsize=12,family=family)
plt.text(4.1,0.322,"\u03C1u=0.98(\u03C1u)$_{US}$",fontsize=12,family=family)
plt.arrow(4.8,0.345,0,-0.01,head_width=0.025, head_length=0.005,fc='k', ec='k')
plt.arrow(4.4,0.330,0,+0.01,head_width=0.025, head_length=0.005,fc='k', ec='k')

plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.savefig('/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp428/Plots/RH428.pdf')
"""

#####################################################
############# Exp 18c ###############################
#####################################################
"""
f="Water"
EoS="HEOS"

P3=4401669.200847257
T3=514.6120306374619
u3=622.0522478762954
Prop3=Cool_PT(P3,T3,f,EoS)
V3=1/Prop3[1]
#cj=fnd_CJS(P3,T3,u3,f,EoS)
Jcj=0.9113570000000001

detCJ=fnd_RHS(P3,T3,u3,Jcj,f,EoS)
Vcj=detCJ[0]
Pcj=detCJ[1]/100000

det=fnd_RHS(P3,T3,u3,1,f,EoS)
Vwd=det[0]
Pwd=det[1]/100000

Vsd=det[3]
Psd=det[4]/100000


rl=RayleighlS(1,P3,T3,u3,f,EoS)
rlcj=RayleighlS(Jcj-0.04,P3,T3,u3,f,EoS)
hlcj=Hugoniotl(P3,T3,0,f,EoS)
P3=P3/100000


plt.figure(figsize=(7.5,5))

family="serif"
axisize=14
lsize=16



font = {'family': family,
        'weight': 'normal',
        'size': 12,
        }

plt.ylabel("Pressure (bar)",fontsize=lsize,family=family)
plt.xlabel("Specific volume (m$^{3}$/kg)",fontsize=lsize,family=family)
#plt.ylim(0.2, 0.4)
#plt.xlim(4, 7)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)



plt.plot(hlcj[0],hlcj[1],"-",color="darkcyan",lw=2,label="Hugoniot")

plt.plot(rl[0],rl[1],"-",color="darkred",lw=1.5,label="Rayleigh")
plt.plot(rlcj[0],rlcj[1],"-",color="darkred",lw=1.5)



plt.plot(V3,P3,"v",ms=8,color="dodgerblue",label="Upstream state")
plt.plot(Vwd,Pwd,"o",ms=7,color="blue",label="Weak detonation")
plt.plot(Vsd,Psd,"o",ms=7,color="darkgreen",label="Strong detonation")
plt.plot(Vcj,Pcj,"^",ms=8,color="black",label="CJ-point")
plt.legend(frameon=False,prop=font)

plt.text(0.010,114,"\u03C1u=(\u03C1u)$_{US}$",fontsize=12,family=family)
plt.text(0.005,88,"\u03C1u=0.91(\u03C1u)$_{US}$",fontsize=12,family=family)
plt.arrow(0.0098,92,0,+3,width=0.0001,head_width=0.0004, head_length=3,fc='k', ec='k')
plt.arrow(0.014,113,0,-3,width=0.0001,head_width=0.0004, head_length=3,fc='k', ec='k')

plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.savefig('/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/18cdprl/Plots18c/RH18c.pdf')
"""

#####################################################
############# Exp 18b ###############################
#####################################################

f="Water"
EoS="HEOS"

P3=3051969.536239191
T3=491.67353717347237
u3=752.8494070533794

Prop3=Cool_PT(P3,T3,f,EoS)
V3=1/Prop3[1]
#cj=fnd_CJS(P3,T3,u3,f,EoS)
Jcj=0.811364768




detCJ=fnd_RHS(P3,T3,u3,Jcj,f,EoS)
Vcj=detCJ[0]
Pcj=detCJ[1]/100000

det=fnd_RHS(P3,T3,u3,1,f,EoS)
Vwd=det[0]
Pwd=det[1]/100000

Vsd=det[3]
Psd=det[4]/100000


rl=RayleighlS(1,P3,T3,u3,f,EoS)
rlcj=RayleighlS1(Jcj-0.075,P3,T3,u3,f,EoS)
hlcj=Hugoniotl(P3,T3,0,f,EoS)
P3=P3/100000


plt.figure(figsize=(7.5,5))

family="serif"
axisize=14
lsize=16



font = {'family': family,
        'weight': 'normal',
        'size': 12,
        }

plt.ylabel("Pressure (bar)",fontsize=lsize,family=family)
plt.xlabel("Specific volume (m$^{3}$/kg)",fontsize=lsize,family=family)
#plt.ylim(0.2, 0.4)
#plt.xlim(4, 7)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)



plt.plot(hlcj[0],hlcj[1],"-",color="darkcyan",lw=2,label="Hugoniot")

plt.plot(rl[0],rl[1],"-",color="darkred",lw=1.5,label="Rayleigh")
plt.plot(rlcj[0],rlcj[1],"-",color="darkred",lw=1.5)



plt.plot(V3,P3,"v",ms=8,color="darkslategrey",label="Upstream state")
plt.plot(Vwd,Pwd,"o",ms=7,color="blue",label="Weak detonation")
plt.plot(Vsd,Psd,"o",ms=7,color="darkgreen",label="Strong detonation")
plt.plot(Vcj,Pcj,"^",ms=8,color="black",label="CJ-point")
plt.legend(frameon=False,prop=font)

plt.text(0.03,81,"\u03C1u=(\u03C1u)$_{US}$",fontsize=12,family=family)
plt.text(0.026,48,"\u03C1u=0.81(\u03C1u)$_{US}$",fontsize=12,family=family)
plt.arrow(0.029,51,0,+3,width=0.0001,head_width=0.0004, head_length=3,fc='k', ec='k')
plt.arrow(0.033,79,0,-3,width=0.0001,head_width=0.0004, head_length=3,fc='k', ec='k')

plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.savefig('/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18b/Plots18b/RH18b.pdf')









