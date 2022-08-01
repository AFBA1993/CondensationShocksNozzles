#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 17:10:18 2020

@author: andres
"""
import matplotlib.pyplot as plt
import pandas as pd



mesh100=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Bierco2results/bierco2100/bierdrpl100_L.csv")
mesh200=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Bierco2results/bierco2200/bierdrpl200_L.csv")
mesh300=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Bierco2results/bierco2300/bierdrpl300_L.csv")
mesh400=pd.read_csv("/media/andres777/DISPOSITIVO//Mestrado/Andressss/Final/Bierco2results/bierco2400/bierdrpl400_L.csv")
mesh500=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Bierco2results/bierco2500/bierdrpl500_L.csv")
mesh1000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Bierco2results/bierco21000/bierdrpl1000_L.csv")
mesh2500=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Bierco2results/bierco22500/bierdrpl2500_L.csv")
mesh5000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Bierco2results/bierco25000/bierdrpl5000_L.csv")
mesh10000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Bierco2results/bierco210000/bierdrpl10000_L.csv")

plt.figure(1)
plt.figure(figsize=(7.5,5))
family="serif"
axisize=14
lsize=16


font = {'family': family,
        'weight': 'normal',
        'size': 12,
        }

font2 = {'family': family,
        'weight': 'normal',
        'size': 13,
        }

plt.xlim(1.2,2)
plt.ylim(0.38,0.42)
plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (mm)",fontsize=lsize,family=family)

y="PPo"
plt.plot(mesh100["x"],mesh100[y],"--",lw=1,color="blue",label="Mesh=100")
plt.plot(mesh200["x"],mesh200[y],"-",lw=1,color="black",label="Mesh=200")
plt.plot(mesh300["x"],mesh300[y],"-",lw=1,color="dodgerblue",label="Mesh=300")
#plt.plot(mesh400["x"],mesh400["PPo"],"-",lw=1,color="grey",label="Mesh=400")
#plt.plot(mesh500["x"],mesh500["PPo"],"-",lw=1,color="mediumseagreen",label="Mesh=500")
plt.plot(mesh1000["x"],mesh1000[y],"--",lw=1,color="forestgreen",label="Mesh=1000")
plt.plot(mesh2500["x"],mesh2500[y],"-",lw=1,color="darkcyan",label="Mesh=2500")
plt.plot(mesh5000["x"],mesh5000[y],"-",lw=1,color="firebrick",label="Mesh=5000")
plt.plot(mesh10000["x"],mesh10000[y],"-",lw=1,color="saddlebrown",label="Mesh=10000")

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)


#plt.grid(color="gainsboro", linestyle='--', linewidth=.75)


plt.legend(frameon=False,prop=font2)
#plt.savefig("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/Bierco2results/Plotsbierco2/meshppocomp.pdf")




"""
drpl=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/Bierco2results/bierco25000/bierdrpl5000_L.csv")
schk=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/Bierco2results/bierco2schk/co2shckbierco2schk_schk.csv")
exp=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Bierco2/co2res/bierexp.csv")
"""
"""
plt.figure(1)
plt.figure(figsize=(7.5,5))
family="serif"
axisize=14
lsize=16

font = {'family': family,
        'weight': 'normal',
        'size': 12,
        }

font2 = {'family': family,
        'weight': 'normal',
        'size': 13,
        }


plt.xlim(0,4)
plt.ylim(0.2,0.6)

plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (mm)",fontsize=lsize,family=family)

plt.plot(drpl["x"],drpl["PPo"],"-",lw=2,color="forestgreen",label="Droplet growth")
plt.plot(schk["x"],schk["PPo"],"-",lw=2,color="darkcyan",label="Weak detonation")
plt.plot(exp["x"],exp["ppo"],"--",lw=2.5,color="black",label="Bier et al. (1990)")

plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)
plt.legend(frameon=False,prop=font2)
plt.savefig("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/Bierco2results/Plotsbierco2/PPoco2comp.pdf")
"""
"""
plt.figure(1)
plt.figure(figsize=(7.5,5))
family="serif"
axisize=14
lsize=16

font = {'family': family,
        'weight': 'normal',
        'size': 12,
        }

font2 = {'family': family,
        'weight': 'normal',
        'size': 13,
        }


plt.xlim(0,2)
plt.ylim(0.0,2E25)

plt.ylabel("Nucleation rate J (m$^{-3}$s$^{-1}$)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (mm)",fontsize=lsize,family=family)

plt.plot(drpl["x"],drpl["J"],"-",lw=2,color="forestgreen",label="Droplet growth")
#plt.plot(schk["x"],schk["PPo"],"-",lw=2,color="darkcyan",label="Weak detonation")
#plt.plot(exp["x"],exp["ppo"],"--",lw=2.5,color="black",label="Bier et al. (1990)")

plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)
plt.legend(frameon=False,prop=font2)
#plt.savefig("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/Bierco2results/Plotsbierco2/Nrateco2.pdf")
"""
"""
plt.figure(2)
plt.figure(figsize=(7.5,5))
family="serif"
axisize=14
lsize=16

font = {'family': family,
        'weight': 'normal',
        'size': 12,
        }

font2 = {'family': family,
        'weight': 'normal',
        'size': 13,
        }


plt.xlim(1,6)
plt.ylim(1e-08,5e-08)

plt.ylabel("Droplet mean radius (m)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (mm)",fontsize=lsize,family=family)

plt.plot(drpl["x"],drpl["rm"],"-",lw=2,color="forestgreen",label="Droplet growth")
#plt.plot(schk["x"],schk["PPo"],"-",lw=2,color="darkcyan",label="Weak detonation")
#plt.plot(exp["x"],exp["ppo"],"--",lw=2.5,color="black",label="Bier et al. (1990)")

plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)
plt.legend(frameon=False,prop=font2)
plt.savefig("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/Bierco2results/Plotsbierco2/rmco2.pdf")
"""