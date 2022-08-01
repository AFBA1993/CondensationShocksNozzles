# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 22:09:22 2019

@author: USER
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def sol_profiles(exp_sphase,exp_2phase,To,Po):
    
    exp=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/lowpress_drpl/Mosesres/Moses417/Mosesexp417.csv")
    sk=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/final2/IPRSV_RESULTS/Moses_iprsv/moseschk/iprsvmoses417_schk.csv")

    
    
    ne_sphase=exp_sphase[15]
    ne_2phase=exp_2phase[0]
    
    c_sphase=0
    c_2phase=0
    
    for i in range (ne_sphase):
        if exp_sphase[4][i]!=0:
            c_sphase=c_sphase+1
        
    for i in range (ne_2phase):
        if exp_2phase[2][i]!=0:
            c_2phase=c_2phase+1
            
    ne=c_sphase+c_2phase
    x_sol=np.zeros(ne)
    P_sol=np.zeros(ne)
    T_sol=np.zeros(ne)
    y_sol=np.zeros(ne)
    S_sol=np.zeros(ne)
    rcrit=np.zeros(ne)
    SubC=np.zeros(ne)
    ei=np.zeros(ne)
    J_sol=np.zeros(ne)
    A_sol=np.zeros(ne)
    Trat=np.zeros(ne)
    Prat=np.zeros(ne)
    u_sol=np.zeros(ne)
    rm_sol=np.zeros(ne)
    rnew_sol=np.zeros(ne)
    emass_sol=np.zeros(ne)
    emom_sol=np.zeros(ne)
    s_sol=np.zeros(ne)
    M_sol=np.zeros(ne)
    eener_sol=np.zeros(ne)
    rgrp1=np.zeros(ne)
    rgrp2=np.zeros(ne)
    delr=np.zeros(ne)
    Nsum=np.zeros(ne)
    Kn=np.zeros(ne)
    Tl=np.zeros(ne)
    Tlnew=np.zeros(ne)
    
    j=0
    
    for i in range (ne):
        if i<=(c_sphase-1):
            x_sol[i]=exp_sphase[0][i]
            P_sol[i]=exp_sphase[2][i]
            T_sol[i]=exp_sphase[4][i]
            S_sol[i]=exp_sphase[9][i]
            J_sol[i]=exp_sphase[10][i]
            A_sol[i]=exp_sphase[16][i]
            u_sol[i]=exp_sphase[7][i]
            M_sol[i]=exp_sphase[8][i]
            s_sol[i]=exp_sphase[6][i]
            SubC[i]=exp_sphase[19][i]
            emass_sol[i]=exp_sphase[22][i]
            emom_sol[i]=exp_sphase[21][i]
            eener_sol[i]=exp_sphase[23][i]
            ei[i]=exp_sphase[24][i]
            #rcrit[i]=exp_sphase[20][i]
            
            j=j+1
        if i>(c_sphase-1):
            x_sol[i]=exp_2phase[1][i-j]
            P_sol[i]=exp_2phase[2][i-j]
            T_sol[i]=exp_2phase[3][i-j]
            y_sol[i]=exp_2phase[5][i-j]
            S_sol[i]=exp_2phase[7][i-j]
            u_sol[i]=exp_2phase[4][i-j]
            J_sol[i]=exp_2phase[8][i-j]
            s_sol[i]=exp_2phase[6][i-j]
            SubC[i]=exp_2phase[16][i-j]
            rm_sol[i]=exp_2phase[9][i-j]
            rnew_sol[i]=exp_2phase[10][i-j]
            emass_sol[i]=exp_2phase[11][i-j]
            emom_sol[i]=exp_2phase[12][i-j]
            ei[i]=exp_2phase[25][i-j]
            M_sol[i]=exp_2phase[26][i-j]
            Tl[i]=exp_2phase[27][i-j]
            Tlnew[i]=exp_2phase[28][i-j]
            eener_sol[i]=exp_2phase[24][i-j]
            rgrp1[i]=exp_2phase[14][i-j]
            rgrp2[i]=exp_2phase[15][i-j]
            rcrit[i]=exp_2phase[17][i-j]
            delr[i]=exp_2phase[18][i-j]
            Nsum[i]=exp_2phase[19][i-j]
            
            Kn[i]=exp_2phase[20][i-j]
       
        Trat[i]=T_sol[i]/To
        Prat[i]=P_sol[i]/Po
        
        
    c_rnew=0
    c_rgrp1=0
    c_rgrp2=0
    c_Tl=0
    
    for j in range (ne):
        if rnew_sol[j]!=0:
            if c_rnew==0:
                xi_rnew=j
            c_rnew+=1
            
    for j in range (ne):
        if rgrp1[j]!=0:
            if c_rgrp1==0:
                xi_rgrp1=j
            c_rgrp1+=1
    
    for j in range (ne):
        if rgrp2[j]!=0:
            if c_rgrp2==0:
                xi_rgrp2=j
            c_rgrp2+=1
            
    for j in range (ne):
        if Tlnew[j]!=0:
            if c_Tl==0:
                xi_Tl=j
            c_Tl+=1
      
    x_rnew=np.zeros(c_rnew)
    y_rnew=np.zeros(c_rnew)
    x_rgrp1=np.zeros(c_rgrp1)
    y_rgrp1=np.zeros(c_rgrp1)
    x_rgrp2=np.zeros(c_rgrp2)
    y_rgrp2=np.zeros(c_rgrp2)
    x_Td=np.zeros(c_Tl)
    y_Td=np.zeros(c_Tl)
    
    for i in range(c_rnew):
        x_rnew[i]=x_sol[xi_rnew]
        y_rnew[i]=rnew_sol[xi_rnew]
        xi_rnew+=1
    
    for i in range(c_rgrp1):
        x_rgrp1[i]=x_sol[xi_rgrp1]
        y_rgrp1[i]=rgrp1[xi_rgrp1]
        xi_rgrp1+=1
    
    for i in range(c_rgrp2):
        x_rgrp2[i]=x_sol[xi_rgrp2]
        y_rgrp2[i]=rgrp2[xi_rgrp2]
        xi_rgrp2+=1
        
    for i in range(c_Tl):
        x_Td[i]=x_sol[xi_Tl]
        y_Td[i]=Tl[xi_Tl]
        xi_Tl+=1
    
    xW=exp_2phase[13]
    PW=exp_2phase[21]
    TW=exp_2phase[22]
    YW=exp_2phase[23]
    
    grid=ne
    
    plt.figure(1)
    plt.plot(x_Td,y_Td)
    plt.plot(x_sol,T_sol)
    plt.plot(sk["x"],sk["T"])
    
    plt.figure(2)
    plt.xlim(8.22,14)
    plt.ylim(0.2,0.6)
    plt.plot(x_sol,Prat)
    plt.plot(exp["x"][:],exp["ppo"][:],"--")
    
    plt.figure(3)
    plt.plot(x_rgrp1,y_rgrp1,"--")
    plt.plot(x_rgrp2,y_rgrp2,"--")
    plt.plot(x_rnew,y_rnew)
    plt.plot(x_sol,rm_sol,"--")
    
            

    """
    plt.figure(1)
    plt.plot(x_sol,emom_sol)
    
    plt.figure(2)
    plt.plot(x_sol,emass_sol)
    
    plt.figure(3)
    plt.plot(x_sol,eener_sol)
    
    plt.figure(4)
    plt.plot(x_sol,s_sol)
    
    plt.figure(5)
    plt.xlim(8.22,14)
    plt.ylim(0.2,0.6)
    plt.plot(x_sol,Prat)
    plt.plot(exp["x"][:],exp["ppo"][:],"--")
    
    plt.figure(6)
    plt.plot(x_sol,S_sol) 
    
    plt.figure(7)
    plt.plot(x_rnew,y_rnew) 
    #plt.plot(expr["x"]/100,expr['rm'],"o")
    
    plt.figure(8)
    plt.plot(x_sol,J_sol) 
    plt.grid()
    
    plt.figure(9)
    plt.plot(x_sol,delr) 
    
    plt.figure(10)
    plt.plot(x_sol,y_sol) 
    
    plt.figure(11)
    plt.plot(x_rgrp1,y_rgrp1)
    
    plt.figure(12)
    plt.plot(x_rgrp2,y_rgrp2)
    
    
    fig, ax1 = plt.subplots()

    ax2 = ax1.twinx()
    ax1.plot(x_sol, Prat, 'g-')
    ax1.plot(exp["x"][:]/10, exp["ppo"], 'o')
    ax1.set_xlim(0.2,0.6)
    ax1.set_ylim(0.25,0.6)
    ax2.plot(x_rnew, y_rnew, 'g-')
    ax2.set_ylim(1e-9,1e-7)
    ax2.set_yscale("log")
    ax2.plot(expr["x"]/10, expr['rm'], 'o')
   """


    plt.show()

      
    
    
    
    return x_sol,P_sol,T_sol,y_sol,S_sol,J_sol,u_sol,rm_sol,rnew_sol,emass_sol,emom_sol,Trat,Prat,xW,PW,TW,YW,SubC,rcrit,delr,Nsum,Kn,x_rgrp1,y_rgrp1,x_rgrp2,y_rgrp2,x_rnew,y_rnew,grid,eener_sol,s_sol,ei,M_sol,Tl,x_Td,y_Td   

def Deliver_res(sol,Folder2send,Name_Data):
    
    resultsL={"x":sol[0],
             "P":sol[1],
             "T":sol[2],
             "u":sol[6],
             "y":sol[3],
             "s":sol[30],
             "rm":sol[7],
             "rnew":sol[8],
             "TTo":sol[11],
             "PPo":sol[12],
             "Emass":sol[9],
             "Emom":sol[10],
             "Eener":sol[29],
             "DelT":sol[17],
             "S":sol[4],
             "J":sol[5],
             "rcrit":sol[18],
             "delr":sol[19],
             "N":sol[20],
             "Kn":sol[21],
             "ei":sol[31],
             "M":sol[32],
             "Tl":sol[33],
           }
    
    
    resultsPoints={"xW":sol[13],"PW":sol[14],"TW":sol[15],"yW":sol[16],"grid":sol[28]}
    resultsg1={"xg1":sol[22],"yg1":sol[23]}
    resultsg2={"xg2":sol[24],"yg2":sol[25]}
    resultsrm={"xnew":sol[26],"ynew":sol[27]}
    resTd={"xTd":sol[34],"yTd":sol[35]}
             
    
    df1=pd.DataFrame(resultsL,columns=["x","P","T","u","y","s","rm","rnew","TTo","PPo","Emass","Emom","Eener","DelT","S",
             "J","rcrit","delr","N","Kn","ei","M","Tl"])
    df2=pd.DataFrame(resultsg1,columns=["xg1","yg1"])
    df3=pd.DataFrame(resultsg2,columns=["xg2","yg2"])
    df4=pd.DataFrame(resultsrm,columns=["xnew","ynew"])
    df5=pd.DataFrame(resultsPoints,columns=["xW","PW","TW","yW","grid"],index=[0])
    df6=pd.DataFrame(resTd,columns=["xTd","yTd"])
    
    df1_to_deliver=Folder2send+Name_Data+"_L.csv"
    df2_to_deliver=Folder2send+Name_Data+"_g1.csv"
    df3_to_deliver=Folder2send+Name_Data+"_g2.csv"
    df4_to_deliver=Folder2send+Name_Data+"_rm.csv"
    df5_to_deliver=Folder2send+Name_Data+"_Points.csv"
    df6_to_deliver=Folder2send+Name_Data+"_tl.csv"
    
    df1.to_csv(df1_to_deliver)
    df2.to_csv(df2_to_deliver)
    df3.to_csv(df3_to_deliver)
    df4.to_csv(df4_to_deliver)
    df5.to_csv(df5_to_deliver)
    df6.to_csv(df6_to_deliver)
    
    return df1_to_deliver, df2_to_deliver, df3_to_deliver, df4_to_deliver, df5_to_deliver,df6_to_deliver

exp=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/lowpress_drpl/Mosesres/Moses410/Mosesexp410.csv")
sk=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/final2/IPRSV_RESULTS/Moses_iprsv/moseschk/iprsvmoses410_schk.csv")


 
#exp=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/lowpress_drpl/Mosesres/Moses417/Mosesexp417.csv")
#sk=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/final2/IPRSV_RESULTS/Moses_iprsv/moseschk/iprsvmoses417_schk.csv")
   
#exp=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/lowpress_drpl/Mosesres/Moses417/Mosesexp417.csv")
#plt.plot(exp["x"][:],exp["ppo"],"-")
