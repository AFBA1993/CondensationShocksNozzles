# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 11:56:04 2019

@author: USER
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from LibiPRSV import iPRSV_mas_PT, iPRSV_mas_PT_sat
from LibiPRSVsat import fnd_Tsat
from Two_phase import Propmix
from Cndshock import cshock
import time
import threading, multiprocessing
import pandas as pd


def diffRH_EQ(P4,P3,T3,u3,q,f):
    prop3=iPRSV_mas_PT(f,P3,T3)
    rho3=1/prop3[1]
    h3=prop3[3]
    
    T4=fnd_Tsat(P4,f)
    prop4=iPRSV_mas_PT_sat(f,P4,T4)
    rhog=1/prop4[2]
    
    X4=(P3-P4+rho3*u3**2.0)*rhog/(rho3*u3)**2.0
    prop4=Propmix(P4,X4,f)
    h4=prop4[3]
    y=-q+h4-h3-0.5*(1/rho3+X4/rhog)*(P4-P3)
    return y

def fndX_RHEQ(P4,P3,T3,u3,f):
    prop3=iPRSV_mas_PT(f,P3,T3)
    rho3=1/prop3[1]
        
    T4=fnd_Tsat(P4,f)
    prop4=iPRSV_mas_PT_sat(f,P4,T4)
    rhog=1/prop4[2]
    
    X4=(P3-P4+rho3*u3**2.0)*rhog/(rho3*u3)**2.0
  

    return X4

def fndPs_RHeq(P3,T3,u3,q,f):
    
    pi=1*P3
    pf=3*P3
    ne=500
    stp=(pf-pi)/ne
    err=np.zeros(ne)
    
    
    for i in range (ne):
        err[i]=diffRH_EQ(pi,P3,T3,u3,q,f)
        print (i)
        if i==(ne-1):
            Pwd=0
        if i >0:
            if (err[i]*err[i-1])<0:
                gssp=pi-0.5*stp
                Pwd=fsolve(diffRH_EQ,gssp,(P3,T3,u3,q,f))
                error=diffRH_EQ(Pwd,P3,T3,u3,q,f)
                if error>1.8e-07:
                    print("error Pwd not found", error)
                
                break
        
        
        pi=pi+stp
        
    
    for i in range (ne):
        err[i]=diffRH_EQ(pi,P3,T3,u3,q,f)
        if i==(ne-1):
            Psd=0
        if i >0:
            if (err[i]*err[i-1])<0:
                gssp=pi-0.5*stp
                Psd=fsolve(diffRH_EQ,gssp,(P3,T3,u3,q,f))
                error=diffRH_EQ(Psd,P3,T3,u3,q,f)
                if error>1.8e-07:
                    print("error Psd not found", error)
                
                break
        pi=pi+stp
        
    return Pwd, Psd

def CJi(q,stp,P3,T3,u3,f):
    qCJ=0
    ne=10000000 
    for i in range (ne):
        Ps=fndPs_RHeq(P3,T3,u3,q,f)
        Pwd=Ps[0]
        Psd=Ps[1]
        diff=Pwd-Psd
        print (diff)
        if Ps[0]==0:
            qnew=q-stp
            stpn=stp/10
            break
        if abs(diff)<104:
            if abs(diff)>0:
                
                qnew=q
                stpn=stp
                qCJ=q
                break
        
        q=q+stp
        
    return qnew,stpn,qCJ,diff


def fnd_CJ(P3,T3,u3,f):
    q=0
    stp=10000
    
    for i in range(10000000000):
        sol=CJi(q,stp,P3,T3,u3,f)
        if sol[2]!=0:
            qCJ=sol[2]
            break
        
        q=sol[0]
        stp=sol[1]
        print (stp,q,sol[3])
        
    Ps=fndPs_RHeq(P3,T3,u3,qCJ,f)
    Pwd_CJ=Ps[0]
    Psd_CJ=Ps[1]
    Xwd_CJ=fndX_RHEQ(Pwd_CJ,P3,T3,u3,f)
    Xsd_CJ=fndX_RHEQ(Psd_CJ,P3,T3,u3,f)
    wd_CJ=Propmix(Pwd_CJ,Xwd_CJ,f)
    sd_CJ=Propmix(Psd_CJ,Xsd_CJ,f)
    Vwd=1/wd_CJ[1]
    Vsd=1/sd_CJ[1]
    swd=wd_CJ[4]
    ssd=sd_CJ[4]
    diffV=Vsd-Vwd
    
    return Pwd_CJ,Xwd_CJ,Vwd,swd,qCJ,Vsd,Psd_CJ,diffV


def Rayleigh(x,P3,rho3,u3,P4,rhog):
    rho4=rhog/x
    y=P4+(rho3*u3)**2.0/rho4-P3-rho3*u3**2.0
    return y

def Rayleigh_rho(P3,rho3,u3,P4,rhog):
    guessx=0.2
    ne=10000
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



def Rayleighl(P3,T3,u3,f):
    prop3=iPRSV_mas_PT(f,P3,T3)
    rho3=1/prop3[1]
    
    
    Pf=2.5*P3
    ne=150
    stp=(Pf-P3)/ne
    P4=P3+stp
    
    Pi=np.zeros(ne)
    Vi=np.zeros(ne)
    
    for i in range (ne):
        
        T4=fnd_Tsat(P4,f)
        prop4=iPRSV_mas_PT(f,P4,T4)
        rhog=1/prop4[1]
        print (rhog)
        rho4=Rayleigh_rho(P3,rho3,u3,P4,rhog)
        Pi[i]=P4
        Vi[i]=1/rho4
        
        P4=P4+stp
        
          

    return Vi,Pi

def diff_Hugoniot(X4,P3,rho3,h3,P4,f,q):
    prop4=Propmix(P4,X4,f)
    h4=prop4[3]
    rho4=prop4[1]
    y=-q+h4-h3-0.5*(1/rho3+1/rho4)*(P4-P3)
    return y

def fnd_Hugoniot(P4,P3,T3,q,f):    


    prop3=iPRSV_mas_PT(f,P3,T3)
    rho3=1/prop3[1]
    h3=prop3[3]
    
       
    xi=0.1
    xf=2
    ne=400
    stpx=(xf-xi)/ne
    
    err=np.zeros(ne)
    for i in range (ne):
        err[i]=diff_Hugoniot(xi,P3,rho3,h3,P4,f,q)
        if i>0:
            if (err[i]*err[i-1])<0:
                xi=xi-0.5*stpx
                X4=fsolve(diff_Hugoniot,xi,(P3,rho3,h3,P4,f,q))
                error=diff_Hugoniot(X4,P3,rho3,h3,P4,f,q)
                if error>1.8e-07:
                    print ("hugoniot no found for",P4,error)
                break
        xi=xi+stpx
    
    prop4=Propmix(P4,X4,f)
    rho4=prop4[1]
    return rho4

def Hugoniotl(P3,T3,q,f):

    Pf=3*P3
    ne=50
    stp=(Pf-P3)/ne
    P4=P3+stp
    
    Pi=np.zeros(ne)
    Vi=np.zeros(ne)
    
    for i in range (ne):
        rhoi=fnd_Hugoniot(P4,P3,T3,q,f)
        Vi[i]=1/rhoi
        Pi[i]=P4/100000
        
        print (P4,i)
        P4=P4+stp
        
    fig,rT=plt.subplots() 
    rT.plot(Vi, Pi)
    rT.set(xlabel="Nozzle Axis", ylabel='Pressure (Pa)')
    rT.grid()
    
    
    

    return Vi,Pi


###############################################################################
###Multiprocessing
def diff_RayleighS(guess_X4,P4,P3,rho3,J3_2,f):
    Prop4=Propmix(P4,guess_X4,f)
    rho4=Prop4[1]
    y=-J3_2+(P4-P3)/(1/rho3-1/rho4)
    return y


def process_RayleighS(guess_X4,stp,ne_p,P4,P3,rho3,J3_2,f,zero):
    err=np.zeros(ne_p)
    
    
    for i in range (ne_p):
        err[i]=diff_RayleighS(guess_X4,P4,P3,rho3,J3_2,f)
        if i>0:
            if (err[i]*err[i-1])<0:
                if abs(err[i]-err[i-1])<5000:
                    guess_X4-=0.5*stp
                    zero[0]=fsolve(diff_RayleighS,guess_X4,(P4,P3,rho3,J3_2,f))
                    break
    
                       
        guess_X4+=stp
        
def fnd_RayleighS_P(P4,P3,rho3,J3_2,f):       
    
        guess_Xi=0.1
        guess_Xf=2
        ne=1000
        stp=(guess_Xf-guess_Xi)/ne
        l=guess_Xf-guess_Xi
        
        nps=12
        
        zero=[0.0]*nps
        
        ne_p=int(ne/nps)+4
        
        processes=[]
        
        for i in range (nps):
            if i ==0:
                guess_Xi=guess_Xi
            if i>0:
                guess_Xi+=l/nps
           
            zero[i]=multiprocessing.Array('d',2) 
            process=multiprocessing.Process(target=process_RayleighS,args=(guess_Xi,stp,ne_p,P4,P3,rho3,J3_2,f,zero[i]))
            process.start()
            processes.append(process)
        
        for process in processes:
                    process.join()
        
       
        res_zeros=np.zeros(nps)
        
        for i in range (nps):
            
            res_zeros[i]=(zero[i][0])
    
        X4=max(res_zeros)
        
        Prop4=Propmix(P4,X4,f)
        
        rho4=Prop4[1]
      
            
        
        return X4,rho4
    
####################################################################
#MUltuprocessing with Thread 

def Rayleigh_threadf(guess_X4t,stp,ne_t,P4,P3,rho3,J3_2,f):
   
    
    global zero, counter_zeros
    counter_zeros=0
    counter_zeros=0
    zero=0.0
    err=[0.0]*ne_t
    for i in range (ne_t):
        if counter_zeros==1:
            break
        err[i]=diff_RayleighS(guess_X4t,P4,P3,rho3,J3_2,f)
        
        if i>0:
            if (err[i]*err[i-1])<0:
                if abs(err[i]-err[i-1])<5000:
                    guess_X4t-=0.5*stp
                    zero=fsolve(diff_RayleighS,guess_X4t,(P4,P3,rho3,J3_2,f))
                    counter_zeros+=1
                    break
    
                       
        guess_X4t+=stp  

    


def process_RayleighS_thread(guess_X4,stp,ne_p,P4,P3,rho3,J3_2,f,root):
    
    nt=1 #number of threads
    l=guess_X4+ne_p*stp
    
    
    ne_t=int(ne_p/nt)+2
    
    threads=[]
    
    for i in range (nt):
        if i ==0:
            guess_X4t=guess_X4
            
        if i>0:
            guess_X4t+=l/nt    
            
        thrd=threading.Thread(target=Rayleigh_threadf,args=(guess_X4t,stp,ne_t,P4,P3,rho3,J3_2,f))
        thrd.start()
        threads.append(thrd)
    
    for thrd in threads:
        thrd.join()
    
    root[0]=zero
        
        
def fnd_RayleighS_main(P4,P3,rho3,J3_2,f):       
    
        guess_Xi=0.001
        guess_Xf=4
        ne=1000
        stp=(guess_Xf-guess_Xi)/ne
        l=guess_Xf-guess_Xi
        
        nps=24 #number of processes
        
        root=[0.0]*nps
        
        ne_p=int(ne/nps)+2
        
        processes=[]
        
        for i in range (nps):
            if i ==0:
                guess_Xi=guess_Xi
            if i>0:
                guess_Xi+=l/nps
           
            root[i]=multiprocessing.Array('d',2) 
            process=multiprocessing.Process(target=process_RayleighS_thread,args=(guess_Xi,stp,ne_p,P4,P3,rho3,J3_2,f,root[i]))
            process.start()
            processes.append(process)
        
        for process in processes:
                    process.join()
        
       
        res_zeros=np.zeros(nps)
        
        for i in range (nps):
            
            res_zeros[i]=(root[i][0])
    
        X4=max(res_zeros)
       
        Prop4=Propmix(P4,X4,f)
        rho4=Prop4[1]
      
            
        
        return X4,rho4    
    

def Rayleigh_linesS(P3,T3,u3,qR,f):
    
    
    Prop3=iPRSV_mas_PT(f,P3,T3)
    rho3=1/Prop3[1]

    J3_2=(rho3*u3)**2.0
    J3_2=qR*J3_2

    Pi=1.03*P3
    Pf=2.5*P3
    ne=50
    
    stp=(Pf-Pi)/ne
    
    P=np.zeros(ne)
    V=np.zeros(ne)
    X=np.zeros(ne)
    for i in range (ne):
        print(i)
        sol=fnd_RayleighS_main(Pi,P3,rho3,J3_2,f)
        X[i]=sol[0]
        V[i]=1/sol[1]
        P[i]=Pi/100000
        Pi+=stp
        
    
    fig,rT=plt.subplots() 
    rT.plot(V, P)
    rT.set(xlabel="Nozzle Axis", ylabel='Pressure (Pa)')
    rT.grid()
    
    return V,P,X
       
###############################################################
#Hugoniot Equation with threads and multiprocessing
    
def fnd_Hugoniot_Thread(Xi,stp,neT,P3,rho3,h3,P4,f,q,T):
    err=[0.0]*neT
    global zero, c
    zero=0.0
    c=0
    for i in range (neT):
        if c>0:
            break
        
        err[i]=diff_Hugoniot(Xi,P3,rho3,h3,P4,f,q)
        if (err[i]*err[i-1])<0:
            Xi-=0.5*stp
            zero=fsolve(diff_Hugoniot,Xi,(P3,rho3,h3,P4,f,q))
            e=diff_Hugoniot(zero,P3,rho3,h3,P4,f,q)
            if abs (e)>1.6e-06:
                print ("Hugoniot not found",e)
            c+=1
            break
        
        Xi+=stp
            
def Hugoniot_thrd(Xi,l,stp,ne_P,P3,rho3,h3,P4,f,q,root):

    NTS=1
    neT=int(ne_P/NTS)+2
    
    threads=[]
    for i in range(NTS):
        if i==0:
            Xi=Xi
        if i>0:
            Xi+=l/NTS
        
        thrd=threading.Thread(target=fnd_Hugoniot_Thread,args=(Xi,stp,neT,P3,rho3,h3,P4,f,q,i))
        thrd.start()
        threads.append(thrd)
        
    for thrd in threads:
        thrd.join()
        
    root[0]=zero
    



def fnd_hugoniot_MP(P4,P3,T3,q,f): 
    
    prop3=iPRSV_mas_PT(f,P3,T3)
    rho3=1/prop3[1]
    h3=prop3[3]
    
    Xi=0.1
    Xf=2
    ne=400
    l=(Xf-Xi)
    stp=(Xf-Xi)/ne
        
    NPS=12 #NUMBRE OF PROCESSES
    ne_P=int(ne/NPS)+2
    
    
    root=[0.0]*NPS
    lT=l/NPS
    processes=[]
    for i in range (NPS):
        if i==0:
            Xi=Xi
        if i>0:
            Xi+=l/NPS
        
        root[i]=multiprocessing.Array("d",1)
        process=multiprocessing.Process(target=Hugoniot_thrd, args=(Xi,lT,stp,ne_P,P3,rho3,h3,P4,f,q,root[i]))
        process.start()
        processes.append(process)
        
    for process in processes:
        process.join()
        
    res_zeros=np.zeros(NPS)
            
    for i in range (NPS):
        res_zeros[i]=(root[i][0])
        
    
        
    X4=max(res_zeros)
    print(diff_Hugoniot(Xi,P3,rho3,h3,P4,f,q))
    Prop4=Propmix(P4,X4,f)
    rho4=Prop4[1]
    
    return rho4,X4

def hugoniotl_MP(P3,T3,q,f):

    Pf=3*P3
    P=1.02*P3
    ne=50
    stp=(Pf-P3)/ne
   
    
    Pi=np.zeros(ne)
    Vi=np.zeros(ne)
    
    for i in range (ne):
        rhoi=fnd_hugoniot_MP(P,P3,T3,q,f)
        Vi[i]=1/rhoi[1]
        Pi[i]=P/100000
        
        #print (P,i)
        P+=stp
        
    
    fig,rT=plt.subplots() 
    rT.plot(Vi, Pi)
    rT.set(xlabel="Nozzle Axis", ylabel='Pressure (Pa)')
    rT.grid()
    

    return Vi,Pi



 
######################################################################## 
########################################################################
#Hugoniot-Rayleigh Multiprocessing
    
##root function
    
def fndX_RHS(guessP4,P3,T3,u3,qR,f):
    prop3=iPRSV_mas_PT(f,P3,T3)
    rho3=1/prop3[1]
    J3=qR*rho3*u3
    J_32=J3*J3
    
    P4=guessP4
    T4=fnd_Tsat(P4,f)
    rhog=iPRSV_mas_PT(f,P4,T4)
    rhog=1/rhog[1]
    X4=((P3-P4)/(qR*J_32)+1/rho3)*rhog
    prop4=Propmix(P4,X4,f)
    V4=1/prop4[1]
    
    return X4 ,V4

def diff_RHS(guessP4,P3,T3,u3,qR,f):
    prop3=iPRSV_mas_PT(f,P3,T3)
    rho3=1/prop3[1]
    h3=prop3[3]
    J3=qR*rho3*u3
    J_32=J3*J3
    
    P4=guessP4
    T4=fnd_Tsat(P4,f)
    rhog=iPRSV_mas_PT(f,P4,T4)
    rhog=1/rhog[1]
    X4=((P3-P4)/(qR*J_32)+1/rho3)*rhog
    prop4=Propmix(P4,X4,f)
    h4=prop4[3]
    y=h4-h3-0.5*(1/rho3+X4/rhog)*(P4-P3)
    return y
    
#Process function
def RHS_P(guessP4,ne_P,stp,P3,T3,u3,qR,f,root):
    c=0
    err=[0.0]*ne_P
    for i in range (ne_P):
        
        err[i]=diff_RHS(guessP4,P3,T3,u3,qR,f)
        if err[i]*err[i-1]<0:
            root[c]=fsolve(diff_RHS,guessP4,(P3,T3,u3,qR,f))
            c+=1
            if c==2:
                break
        guessP4+=stp

#Process main function
def fnd_RH_MP(P3,T3,u3,qR,f):
    
    Pi=1.01*P3
    Pf=3*P3
    ne=600
    stp=(Pf-Pi)/ne
    l=Pf-Pi
    NPS=24
    ne_P=int(ne/NPS)+4
    
    
    
    process=[]
    root=[0.0]*NPS
    for i in range (NPS):
        if i==0:
            Pi=Pi
            
        if i>0:
            Pi+=l/NPS
            
        Pf=Pi+stp*ne_P
        
        root[i]=multiprocessing.Array("d",2)
        p=multiprocessing.Process(target=RHS_P,args=(Pi,ne_P,stp,P3,T3,u3,qR,f,root[i]))
        p.start()
        process.append(p)
    
    for p in process:
        p.join()
        
     
    res_zeros=[0.0]*NPS
    for i in range (NPS):
        res_zeros[i]=root[i][:]
    
    
    Psd=max(res_zeros)
    Psd=max(Psd)
    if Psd>0:
    
        ePsd=diff_RHS(Psd,P3,T3,u3,qR,f)
        if abs(ePsd)>1.8e-08:
            print("Psd not found",ePsd)
        
        for i in range(NPS):
            if res_zeros[i][0]==0:
                res_zeros[i][0]=Psd
            if res_zeros [i][1]==0:
                res_zeros [i][1]=Psd
        Pwd=min(res_zeros)
        Pwd=min(Pwd)
        ePwd=diff_RHS(Pwd,P3,T3,u3,qR,f)
        if abs(ePwd)>1.8e-08:
            print("Pwd not found",ePwd) 
        
        propsd=fndX_RHS(Psd,P3,T3,u3,qR,f)
        Xsd=propsd[0]
        Vsd=propsd[1]
        
        propwd=fndX_RHS(Pwd,P3,T3,u3,qR,f)
        Xwd=propwd[0]
        Vwd=propwd[1]
        
    if Psd==0:
        Pwd=0
        Psd=0
        Xwd=0
        Xsd=0
        Vsd=0
        Vwd=0
    return Pwd, Xwd, Vwd,  Psd, Xsd, Vsd







def CJS(gssPR,stpi,P3,T3,u3,f):
    dd=np.zeros(10000000)
    lol1=fnd_RH_MP(P3,T3,u3,gssPR,f)
    tol=0.1
    j=0
    for i in range (1000):
        lol1=fnd_RH_MP(P3,T3,u3,gssPR,f)
        Pwd=lol1[0]
        Vwd=lol1[2]
        Xwd=lol1[1]
        
        Psd=lol1[3]
        Vsd=lol1[5]
        Xsd=lol1[4]
      
        d=abs(Pwd-Psd)
        dd[i]=d
        if d==0:
            gssPR=gssPR-stpi
            stpi=stpi/10
            VCJ=0
            PCJ=0
            qR=0
            break
        if i>0:
            
            l=abs(dd[i]-dd[i-1])
            if (l)<tol:
               if l>0.0:
                    VCJ=Vwd
                    PCJ=lol1[0]
                    qR=gssPR
                    j=1
                    break
            print (d,gssPR,l)
        if d<7:
            if d>0:
                VCJ=Vwd
                PCJ=lol1[0]
                qR=gssPR
                break
        gssPR=gssPR+stpi
    return gssPR, stpi, VCJ, PCJ, d, qR, j

def fnd_CJS(P3,T3, u3, f):
    
    gssPR=1
    stpi=-0.1
    
    
    for i in range (20):
        cj=CJS(gssPR,stpi,P3,T3,u3,f)
        gssPR=cj[0]
        stpi=cj[1]
        
        d=cj[4]
        if cj[2]!=0:
            VCJ=cj[2]
            PCJ=cj[3]
            JPR=cj[5]
            break
        print(i)
        
        
        if cj[6]==1:
                        
            VCJ=cj[2]
            PCJ=cj[3]
            JPR=cj[5]
           
            break
        
    return VCJ, PCJ, JPR, d

"""
f="Water"
    
P3=27159.052712507626
T3=299.3732780073999
u3=577.0779191826138
#qR=0.9761

t1=time.time()
c=fnd_CJS(P3,T3, u3, f)
Vcj=c[0]
Pcj=c[1]/100000
qR=c[2]

rlc=Rayleigh_linesS(P3,T3,u3,1,f)

ck=fnd_RH_MP(P3,T3,u3,1,f)
Vwd=ck[2]
Pwd=ck[0]/100000

ck2=cshock(P3,T3,u3,0,f)
Vwd1=1/ck2[6]
Pwd1=ck2[5]

Vsd1=1/ck2[9]
Psd1=ck2[10]

Vsd=ck[5]
Psd=ck[3]/100000
rl=Rayleigh_linesS(P3,T3,u3,qR-0.04,f)
hl=Hugoniotl(P3,T3,0,f)


plt.plot(Vcj,Pcj,'o')
plt.plot(Vsd,Psd,'o')
plt.plot(Vwd,Pwd,'o')
plt.plot(rl[0],rl[1],'-')
plt.plot(rlc[0],rlc[1],'-')
plt.plot(hl[0],hl[1],'-')

t2=time.time()-t1

"""

#############################################################
#IPRSV1074###################################################


"""
f="Water"
    
P3=30893.414004972496
T3=302.73764277904104
u3=539.33861317976812


#c=fnd_CJS(P3,T3, u3, f)

#rl=Rayleigh_linesS(P3,T3,u3,1,f)
#hlcj=Hugoniotl(P3,T3,0,f)

#rhlines={"rv":rl[0],"rp":rl[1],"hv":hlcj[0],"hp":hlcj[1]}
#df=pd.DataFrame(rhlines,columns=["rv","rp","hv","hp"])
#df.to_csv("/media/multiusers/DISPOSITIVO/mestrado/Andressss/JHeinze2015/iPRSVcond/iprsv1074RH.csv")

rhlines=pd.read_csv("/media/multiusers/DISPOSITIVO/mestrado/Andressss/JHeinze2015/iPRSVcond/iprsv1074RH.csv")
hv=rhlines["hv"]
hp=rhlines["hp"]

rv=rhlines["rv"]
rp=rhlines["rp"]

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
plt.plot(hv,hp,"-",color="darkcyan",lw=2,label="Hugoniot")
plt.plot(rv,rp,"-",color="darkred",lw=1.5,label="Rayleigh")
plt.legend(frameon=False,prop=font)
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.savefig('RHiprsv1074.pdf')
"""

################################################
#####IPRSV 112.2
"""

f="Water"
    
P3=28595.784555330261
T3=300.82500104807178
u3=562.333294387514

V3=iPRSV_mas_PT(f,P3,T3)
V3=V3[1]
#c=fnd_CJS(P3,T3, u3, f)
Vcj=3.76903439
Pcj=42064.7565533707/100000
Jcj=0.9775664790000003

#rl=Rayleigh_linesS(P3,T3,u3,1,f)
#hlcj=Hugoniotl(P3,T3,0,f)

#rhlines={"rv":rl[0],"rp":rl[1],"hv":hlcj[0],"hp":hlcj[1]}
#f=pd.DataFrame(rhlines,columns=["rv","rp","hv","hp"])
#df.to_csv("/media/multiusers/DISPOSITIVO/mestrado/Andressss/JHeinze2015/iPRSVcond/iprsv1122RH.csv")

rhlines=pd.read_csv("/media/multiusers/DISPOSITIVO/mestrado/Andressss/JHeinze2015/iPRSVcond/iprsv1122RH.csv")    


hv=rhlines["hv"]
hp=rhlines["hp"]

rv=rhlines["rv"]
rp=rhlines["rp"]


#rlcj=Rayleigh_linesS(P3,T3,u3,Jcj-0.04,f)
#rlcj={"rv":rlcj[0],"rp":rlcj[1]}
#df=pd.DataFrame(rlcj,columns=["rv","rp"])
df.to_csv("/media/multiusers/DISPOSITIVO/mestrado/Andressss/JHeinze2015/iPRSVcond/iprsv1122RHcj.csv")
rlcj=pd.read_csv("/media/multiusers/DISPOSITIVO/mestrado/Andressss/JHeinze2015/iPRSVcond/iprsv1122RHcj.csv") 


rlcjv=rlcj["rv"]
rlcjp=rlcj["rp"]


det=fnd_RH_MP(P3,T3,u3,1,f)
Pwd=det[0]/100000
Vwd=det[2]

Psd=det[3]/100000
Vsd=det[5]

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
plt.xlim(3, 4.88)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)



plt.plot(hv,hp,"-",color="darkcyan",lw=2,label="Hugoniot")
plt.plot(rv,rp,"-",color="darkred",lw=1.5,label="Rayleigh")
plt.plot(rlcjv,rlcjp,"-",color="darkred",lw=1.5,label="")



plt.plot(V3,P3,"v",ms=8,color="cornflowerblue",label="State 3")
plt.plot(Vwd,Pwd,"D",ms=5,color="forestgreen",label="Weak detonation")
plt.plot(Vsd,Psd,"D",ms=5,color="blue",label="Strong detonation")
plt.plot(Vcj,Pcj,"o",ms=6,color="black",label="CJ-point")
plt.legend(frameon=False,prop=font)
plt.text(3.05,0.44,"J=0.9775J$_{3}$",fontsize=12,family=family)
plt.text(3.55,0.50,"J=J$_{3}$",fontsize=12,family=family)
plt.arrow(3.6,0.498,0,-0.03,head_width=0.03, head_length=0.01,fc='k', ec='k')
plt.arrow(3.225,0.455,0,+0.02,head_width=0.032, head_length=0.012,fc='k', ec='k')
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.savefig('RHiprsv112.pdf')
"""

#####################################################
#####115.2
"""

f="Water"
    
P3=27251.618212849382
T3=299.62565580590365
u3=576.26564831171561

V3=iPRSV_mas_PT(f,P3,T3)
V3=V3[1]
#c=fnd_CJS(P3,T3, u3, f)
Vcj=3.93421059
Pcj=40177.92716194912/100000
Jcj=0.9611639268999997


#rl=Rayleigh_linesS(P3,T3,u3,1,f)
#hlcj=Hugoniotl(P3,T3,0,f)

#rhlines={"rv":rl[0],"rp":rl[1],"hv":hlcj[0],"hp":hlcj[1]}
#df=pd.DataFrame(rhlines,columns=["rv","rp","hv","hp"])
#df.to_csv("/media/multiusers/DISPOSITIVO/mestrado/Andressss/JHeinze2015/iPRSVcond/iprsv1152RH.csv")

rhlines=pd.read_csv("/media/multiusers/DISPOSITIVO/mestrado/Andressss/JHeinze2015/iPRSVcond/iprsv1152RH.csv")    

hv=rhlines["hv"]
hp=rhlines["hp"]

rv=rhlines["rv"]
rp=rhlines["rp"]


#rlcj=Rayleigh_linesS(P3,T3,u3,Jcj-0.07,f)
#rlcj={"rv":rlcj[0],"rp":rlcj[1]}
#df=pd.DataFrame(rlcj,columns=["rv","rp"])
#df.to_csv("/media/multiusers/DISPOSITIVO/mestrado/Andressss/JHeinze2015/iPRSVcond/iprsv1152RHcj.csv")
#rlcj=pd.read_csv("/media/multiusers/DISPOSITIVO/mestrado/Andressss/JHeinze2015/iPRSVcond/iprsv1152RHcj.csv") 


rlcjv=rlcj["rv"]
rlcjp=rlcj["rp"]



det=fnd_RH_MP(P3,T3,u3,1,f)
Pwd=det[0]/100000
Vwd=det[2]

Psd=det[3]/100000
Vsd=det[5]

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
plt.ylim(0.26, 0.55)
plt.xlim(3, 5.1)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)



plt.plot(hv,hp,"-",color="darkcyan",lw=2,label="Hugoniot")
plt.plot(rv,rp,"-",color="darkred",lw=1.5,label="Rayleigh")
plt.plot(rlcjv,rlcjp,"-",color="darkred",lw=1.5,label="")



plt.plot(V3,P3,"v",ms=8,color="crimson",label="State 3")
plt.plot(Vwd,Pwd,"D",ms=5,color="forestgreen",label="Weak detonation")
plt.plot(Vsd,Psd,"D",ms=5,color="blue",label="Strong detonation")
plt.plot(Vcj,Pcj,"o",ms=6,color="black",label="CJ-point")
plt.legend(frameon=False,prop=font)
plt.text(3.05,0.43,"J=0.96116J$_{3}$",fontsize=12,family=family)
plt.text(3.55,0.52,"J=J$_{3}$",fontsize=12,family=family)
plt.arrow(3.6,0.505,0,-0.03,head_width=0.03, head_length=0.01,fc='k', ec='k')
plt.arrow(3.225,0.448,0,+0.02,head_width=0.032, head_length=0.012,fc='k', ec='k')
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.savefig('RHiprsv1152.pdf')
"""

#####################################################
#####119.7

"""
f="Water"

P3=25412.242025784853
T3=297.97764661261124
u3=596.10303331846683

V3=iPRSV_mas_PT(f,P3,T3)
V3=V3[1]
#c=fnd_CJS(P3,T3, u3, f)

Vcj=4.18477729
Pcj=37613.97741073121/100000
Jcj=0.9387700000000001



#rl=Rayleigh_linesS(P3,T3,u3,1,f)
#hlcj=Hugoniotl(P3,T3,0,f)

#rhlines={"rv":rl[0],"rp":rl[1],"hv":hlcj[0],"hp":hlcj[1]}
#df=pd.DataFrame(rhlines,columns=["rv","rp","hv","hp"])
#df.to_csv("/media/multiusers/DISPOSITIVO/mestrado/Andressss/JHeinze2015/iPRSVcond/iprsv1197RH.csv")

rhlines=pd.read_csv("/media/multiusers/DISPOSITIVO/mestrado/Andressss/JHeinze2015/iPRSVcond/iprsv1197RH.csv")    

hv=rhlines["hv"]
hp=rhlines["hp"]

rv=rhlines["rv"]
rp=rhlines["rp"]


#rlcj=Rayleigh_linesS(P3,T3,u3,Jcj-0.11,f)
#rlcj={"rv":rlcj[0],"rp":rlcj[1]}
#df=pd.DataFrame(rlcj,columns=["rv","rp"])
#df.to_csv("/media/multiusers/DISPOSITIVO/mestrado/Andressss/JHeinze2015/iPRSVcond/iprsv1197RHcj.csv")
rlcj=pd.read_csv("/media/multiusers/DISPOSITIVO/mestrado/Andressss/JHeinze2015/iPRSVcond/iprsv1197RHcj.csv") 


rlcjv=rlcj["rv"]
rlcjp=rlcj["rp"]




det=fnd_RH_MP(P3,T3,u3,1,f)
Pwd=det[0]/100000
Vwd=det[2]

Psd=det[3]/100000
Vsd=det[5]

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
plt.ylim(0.245, 0.55)
plt.xlim(2.8, 5.5)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)



plt.plot(hv,hp,"-",color="darkcyan",lw=2,label="Hugoniot")
plt.plot(rv,rp,"-",color="darkred",lw=1.5,label="Rayleigh")
plt.plot(rlcjv,rlcjp,"-",color="darkred",lw=1.5,label="")



plt.plot(V3,P3,"v",ms=8,color="dodgerblue",label="State 3")
plt.plot(Vwd,Pwd,"D",ms=5,color="forestgreen",label="Weak detonation")
plt.plot(Vsd,Psd,"D",ms=5,color="blue",label="Strong detonation")
plt.plot(Vcj,Pcj,"o",ms=6,color="black",label="CJ-point")
plt.legend(frameon=False,prop=font)
plt.text(3.05,0.42,"J=0.93877J$_{3}$",fontsize=12,family=family)
plt.text(3.55,0.52,"J=J$_{3}$",fontsize=12,family=family)
plt.arrow(3.6,0.515,0,-0.03,head_width=0.03, head_length=0.01,fc='k', ec='k')
plt.arrow(3.225,0.44,0,+0.02,head_width=0.032, head_length=0.012,fc='k', ec='k')
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.savefig('RHiprsv1197.pdf')
"""






















