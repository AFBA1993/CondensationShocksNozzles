# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 22:09:22 2019

@author: USER
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
def sol_profiles(exp_sphase,exp_2phase,To,Po):
    
    ne_sphase=exp_sphase[15]
    ne_2phase=exp_2phase[12]
    
    c_sphase=0
    c_2phase=0
    
    for i in range (ne_sphase):
        if exp_sphase[4][i]!=0:
            c_sphase=c_sphase+1
        
    for i in range (ne_2phase):
        if exp_2phase[1][i]!=0:
            c_2phase=c_2phase+1
            
    ne=c_sphase+c_2phase
    x_sol=np.zeros(ne)
    P_sol=np.zeros(ne)
    T_sol=np.zeros(ne)
    y_sol=np.zeros(ne)
    S_sol=np.zeros(ne)
    s_sol=np.zeros(ne)
    J_sol=np.zeros(ne)
    V_sol=np.zeros(ne)
    A_sol=np.zeros(ne)
    Trat=np.zeros(ne)
    Prat=np.zeros(ne)
    rho_sol=np.zeros(ne)
    u_sol=np.zeros(ne)
    M_sol=np.zeros(ne)
    e_sol=np.zeros(ne)
    nei=np.zeros(ne)
    ener=np.zeros(ne)
    emom=np.zeros(ne)
    emass=np.zeros(ne)
    mass=np.zeros(ne)
    j=0
    
    for i in range (ne):
        if i<=(c_sphase-1):
            x_sol[i]=exp_sphase[0][i]
            P_sol[i]=exp_sphase[2][i]
            T_sol[i]=exp_sphase[4][i]
            S_sol[i]=exp_sphase[9][i]
            J_sol[i]=exp_sphase[10][i]
            V_sol[i]=1/exp_sphase[3][i]
            A_sol[i]=exp_sphase[16][i]
            rho_sol[i]=exp_sphase[3][i]
            u_sol[i]=exp_sphase[7][i]
            e_sol[i]=exp_sphase[22][i]
            M_sol[i]=exp_sphase[8][i]
            ener[i]=exp_sphase[21][i]
            emom[i]=exp_sphase[19][i]
            emass[i]=exp_sphase[20][i]
            s_sol[i]=exp_sphase[6][i]
            nei[i]=ne
            
            j=j+1
        if i>(c_sphase-1):
            x_sol[i]=exp_2phase[0][i-j]
            P_sol[i]=exp_2phase[1][i-j]
            T_sol[i]=exp_2phase[3][i-j]
            y_sol[i]=exp_2phase[13][i-j]
            V_sol[i]=1/exp_2phase[2][i-j]
            rho_sol[i]=exp_2phase[2][i-j]
            S_sol[i]=1
            e_sol[i]=exp_2phase[9][i-j]
            M_sol[i]=exp_2phase[7][i-j]
            s_sol[i]=exp_2phase[5][i-j]
            A_sol[i]=exp_2phase[14][i-j]
            u_sol[i]=exp_2phase[6][i-j]
            mass[i]=exp_2phase[8][i-j]
            emom[i]=exp_2phase[15][i-j]
            emass[i]=exp_2phase[16][i-j]
            ener[i]=exp_2phase[17][i-j]
        Trat[i]=T_sol[i]/To
        Prat[i]=P_sol[i]/Po
            
       
            
       
    plt.figure(1)
    plt.plot(x_sol,emom)
    
    plt.figure(2)
    plt.plot(x_sol,emass)
    
    plt.figure(3)
    plt.plot(x_sol,ener)
    
    plt.figure(4)
    #plt.xlim(8.22,14)
    #plt.ylim(0.2,0.6)
    plt.plot(x_sol,Prat)
    
            
    """
       
    fig,rT=plt.subplots() 
    rT.plot(x_sol, P_sol)
    rT.set(xlabel="Nozzle Axis", ylabel='Pressure (Pa)')
    rT.grid()
    
    
    fig,rT=plt.subplots() 
    rT.plot(x_sol, T_sol)
    rT.set(xlabel="Nozzle Axis", ylabel='TemperatUre (K)')
    rT.grid()
    
    
    fig,rT=plt.subplots() 
    rT.plot(x_sol, y_sol)
    rT.set(xlabel="Nozzle Axis", ylabel='Wetness mass fraction (-)')
    rT.grid()
    
    
    fig,rT=plt.subplots() 
    rT.plot(x_sol, S_sol)
    rT.set(xlabel="Nozzle Axis", ylabel='Supersaturation (K)')
    rT.grid()
    
    fig,rT=plt.subplots() 
    rT.plot(x_sol, J_sol)
    rT.set(xlabel="Nozzle Axis", ylabel='Nucleation Rate (K)')
    rT.grid()
    
    
    fig,rT=plt.subplots() 
    rT.plot(x_sol, P_sol)
    rT.set(x_label="Nozzle Axis", ylabel='Pressure (Pa)')
    rT.grid()
    """

    return x_sol, P_sol,T_sol,y_sol,S_sol,J_sol,V_sol,A_sol,rho_sol, u_sol,Trat,Prat,s_sol,M_sol,e_sol,emom,ener,emass,mass,nei 

def Deliver_res(sol,Folder2send,Name_Data):  


    results={'x':sol[0],
          "A":sol[7],
          "P":sol[1],
          "rho":sol[8],
          "u":sol[9],
          "y":sol[3],
          "PPo":sol[11],
          "TTo":sol[10],
          "T":sol[2],
          "M":sol[13],
          "e":sol[14],
          "emom":sol[15],
          "eener":sol[16],
          "emass":sol[17],
          "mass":sol[18],
          "grid":sol[19],
          "s":sol[12]}


    df1=pd.DataFrame(results,columns=["x","A","P","rho","u","y","PPo","TTo","T","M","e","emom","eener","emass","mass","grid","s"])
    df1_to_deliver=Folder2send+Name_Data+"_schk.csv"
    df1.to_csv(df1_to_deliver)
    return df1_to_deliver