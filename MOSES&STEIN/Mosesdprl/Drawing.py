#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 17:10:18 2020

@author: andres
"""
import matplotlib.pyplot as plt
import pandas as pd

#####Mesh validation####

"""
mesh100=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/ResultsMosescomp/Mesh410/moses100/410Gyardrp100moses_L.csv")
mesh200=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/ResultsMosescomp/Mesh410/Moses200/410Gyardrp200moses_L.csv")
mesh300=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/ResultsMosescomp/Mesh410/Moses300/410Gyardrp300moses_L.csv")
mesh400=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/ResultsMosescomp/Mesh410/Moses400/410Gyardrp400moses_L.csv")
mesh500=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/ResultsMosescomp/Mesh410/Moses500/410Gyardrp500moses_L.csv")
mesh1000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/ResultsMosescomp/Mesh410/Moses1000/410Gyardrp1000moses_L.csv")
mesh2500=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/ResultsMosescomp/Mesh410/Moses2500/410Gyardrp2500moses_L.csv")
mesh4000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/ResultsMosescomp/Mesh410/Moses4000/410Gyardrp4000moses_L.csv")
mesh5000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/ResultsMosescomp/Mesh410/Moses5000/410Gyardrp5000moses_L.csv")
mesh10000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/ResultsMosescomp/Mesh410/Moses10000/410Gyardrp10000moses_L.csv")





#Pressure profile


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



xi="\u03BE"
psi="\u03C8"

plt.xlim(10,12)
plt.ylim(0.36,0.42)

plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (cm)",fontsize=lsize,family=family)

y="PPo"

plt.plot(mesh100["x"]*100,mesh100[y],"--", color="darkred",lw=1.5,label="Mesh=100")
plt.plot(mesh200["x"]*100,mesh200[y],"--", color="steelblue",lw=1.5,label="Mesh=200")
plt.plot(mesh300["x"]*100,mesh300[y],"--", color="black",lw=1.5,label="Mesh=300")
#plt.plot(mesh400["x"],mesh400["PPo"],"-",lw=1.0,color="black",label="Mesh=2500")
#plt.plot(mesh500["x"]*100,mesh500[y],"-",lw=1.0,color="forestgreen",label="Mesh=500")
plt.plot(mesh1000["x"]*100,mesh1000[y],"-",color="darkcyan",lw=1.2,label="Mesh=1000")
#plt.plot(mesh2000["x"]/10,mesh2000[y],"--",color="green",lw=1.2,label="Mesh=2000")
#plt.plot(mesh3000["x"],mesh3000["PPo"],"-",color="darkgreen",lw=1.2,label="Mesh=3000")
plt.plot(mesh4000["x"]*100,mesh4000[y],"-",color="seagreen",lw=1.2,label="Mesh=4000")
plt.plot(mesh5000["x"]*100,mesh5000[y],"-",color="darkblue",lw=1.0,label="Mesh=5000")
plt.plot(mesh10000["x"]*100,mesh10000[y],"-",color="black",lw=1.0,label="Mesh=10000")

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

#plt.savefig("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/ResultsMosescomp/exp410/Plots/mesh410.pdf")
"""








####################################################################
####################################################################
####################################################################
#EXPERIMENT410

#shock=pd.read_csv("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp410/410shck/410schk_schk.csv")
#drpl=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsMosescomp/exp410/410drpl/410drpl410Gyardrp_L.csv")


#drplr=pd.read_csv("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp410/410drpl/410drpl410Gyardrp_L.csv")
#exp=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/lowpress_drpl/Mosesres/Moses410/Mosesexp410.csv")


##############################################################################
#Pressure profile

"""
plt.figure(1)
plt.figure(figsize=(7,5))
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



xi="\u03BE"
psi="\u03C8"

plt.xlim(8.22/100,14/100)
plt.ylim(0.2,0.7)

plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

#plt.plot(shock["x"]/100,shock["PPo"],"-",lw=2.5,color="forestgreen",label="Weak detonation")
plt.plot(drpl["x"],drpl["PPo"],"-",lw=2.5,color="dodgerblue",label="Droplet growth")
#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
plt.plot(exp["x"]/100,exp["ppo"],"--",lw=2.5,color="black",label="Moses and Stein (1978)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
#plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp410/Plots/PPocomp.pdf")
"""


#Nucleation rate
"""
plt.figure(1)
plt.figure(figsize=(7,5))
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



xi="\u03BE"
psi="\u03C8"

plt.xlim(5/100,14/100)
plt.ylim(1,1e25)
plt.yscale("log")


plt.ylabel("Nucleation rate J (m$^{-3}$s$^{-1}$)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

#plt.plot(shock["x"]/100,shock["PPo"],"-",lw=2.5,color="forestgreen",label="Weak detonation")
plt.plot(drpl["x"],drpl["J"],"-",lw=2.5,color="dodgerblue",label="Droplet growth")
#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
#plt.plot(exp["x"]/100,exp["ppo"],"--",lw=2.5,color="black",label="Moses and Stein (1978)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
#plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp410/Plots/nucl410log.pdf")
"""

#Radius
"""
plt.figure(1)
plt.figure(figsize=(7,5))
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



xi="\u03BE"
psi="\u03C8"

plt.xlim(10/100,16/100)
plt.ylim(0.25e-08,1.75e-08)
#plt.yscale("log")


plt.ylabel("Droplet mean radius (m)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

#plt.plot(shock["x"]/100,shock["PPo"],"-",lw=2.5,color="forestgreen",label="Weak detonation")
plt.plot(drplr["x"],drplr["rm"],"-",lw=2.5,color="dodgerblue",label="Droplet growth")
#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
#plt.plot(exp["x"]/100,exp["ppo"],"--",lw=2.5,color="black",label="Moses and Stein (1978)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

#plt.legend(frameon=False,prop=font2)
plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp410/Plots/rm410.pdf")
"""
###############################################################
###############################################################
###############################################################
#EXPERIMENT417

#shock=pd.read_csv("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp417/417shck/417schk_schk.csv")
#drpl=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsMosescomp/exp417/417drpl/417Gyardrp_L.csv")
#exp=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/lowpress_drpl/Mosesres/Moses417/Mosesexp417.csv")

#Pressure profile

"""
plt.figure(1)
plt.figure(figsize=(7,5))
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



xi="\u03BE"
psi="\u03C8"

plt.xlim(8.22/100,14/100)
plt.ylim(0.2,0.7)

plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

#plt.plot(shock["x"]/100,shock["PPo"],"-",lw=2.5,color="royalblue",label="Weak detonation")
plt.plot(drpl["x"],drpl["PPo"],"-",lw=2.5,color="darkcyan",label="Droplet growth")
plt.plot(exp["x"]/100,exp["ppo"],"--",lw=2.5,color="black",label="Moses and Stein (1978)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
#plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp417/Plots/PPo417.pdf")


"""

"""
#Nucleation rate

plt.figure(1)
plt.figure(figsize=(7,5))
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



xi="\u03BE"
psi="\u03C8"

plt.xlim(8/100,14/100)
plt.ylim(1,8e22)
#plt.yscale("log")


plt.ylabel("Nucleation rate J (m$^{-3}$s$^{-1}$)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

#plt.plot(shock["x"]/100,shock["PPo"],"-",lw=2.5,color="forestgreen",label="Weak detonation")
plt.plot(drpl["x"],drpl["J"],"-",lw=2.5,color="darkcyan",label="Droplet growth")
#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
#plt.plot(exp["x"]/100,exp["ppo"],"--",lw=2.5,color="black",label="Moses and Stein (1978)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)


plt.legend(frameon=False,prop=font2)
plt.savefig("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsMosescomp/exp417/Plots/Nrate417.pdf")
"""

"""
#Radius

plt.figure(1)
plt.figure(figsize=(7,5))
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



xi="\u03BE"
psi="\u03C8"

plt.xlim(10/100,16/100)
plt.ylim(0.35e-08,1.75e-08)
#plt.yscale("log")


plt.ylabel("Droplet mean radius (m)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

#plt.plot(shock["x"]/100,shock["PPo"],"-",lw=2.5,color="forestgreen",label="Weak detonation")
plt.plot(drpl["x"],drpl["rnew"],"-",lw=2.5,color="darkcyan",label="Droplet growth")
#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
#plt.plot(exp["x"]/100,exp["ppo"],"--",lw=2.5,color="black",label="Moses and Stein (1978)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp417/Plots/rm417.pdf")
"""
#plt.legend(frameon=False,prop=font2)

########################################################################################33
##############################################################################################
##############################################################################################
####EXP428
"""
#shock=pd.read_csv("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp428/428schk/428schk_schk.csv")
drpl=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsMosescomp/exp428/428dprl/428Gyardrp_L.csv")
exp=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/lowpress_drpl/Mosesres/Moses428/Mosesexp428.csv")

#PPo

plt.figure(1)
plt.figure(figsize=(7,5))
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



xi="\u03BE"
psi="\u03C8"

plt.xlim(8.22/100,14/100)
plt.ylim(0.2,0.7)

plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

#plt.plot(shock["x"]/100,shock["PPo"],"-",lw=2.5,color="limegreen",label="Weak detonation")
plt.plot(drpl["x"],drpl["PPo"],"-",lw=2.5,color="darkslategrey",label="Droplet growth")
plt.plot(exp["x"]/100,exp["ppo"],"--",lw=2.5,color="black",label="Moses and Stein (1978)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
#plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp428/Plots/PPo428.pdf")

#Nucleation rate
"""
"""
plt.figure(1)
plt.figure(figsize=(7,5))
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


xi="\u03BE"
psi="\u03C8"

plt.xlim(8.22/100,13/100)
plt.ylim(1,0.16e24)
plt.yscale("log")


plt.ylabel("Nucleation rate J (m$^{-3}$s$^{-1}$)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

#plt.plot(shock["x"]/100,shock["PPo"],"-",lw=2.5,color="forestgreen",label="Weak detonation")
plt.plot(drpl["x"],drpl["J"],"-",lw=2.5,color="darkslategrey",label="Droplet growth")
#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
#plt.plot(exp["x"]/100,exp["ppo"],"--",lw=2.5,color="black",label="Moses and Stein (1978)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)


plt.legend(frameon=False,prop=font2)
plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp428/Plots/Nrate428log.pdf")
"""

#Radius
"""
plt.figure(1)
plt.figure(figsize=(7,5))
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



xi="\u03BE"
psi="\u03C8"

plt.xlim(10/100,16/100)
plt.ylim(0.27e-08,1.3e-08)
#plt.yscale("log")


plt.ylabel("Droplet mean radius (m)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

#plt.plot(shock["x"]/100,shock["PPo"],"-",lw=2.5,color="forestgreen",label="Weak detonation")
plt.plot(drpl["x"],drpl["rnew"],"-",lw=2.5,color="darkslategrey",label="Droplet growth")
#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
#plt.plot(exp["x"]/100,exp["ppo"],"--",lw=2.5,color="black",label="Moses and Stein (1978)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)
plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp428/Plots/rm428.pdf")
"""
#######################
###################################
################################
#EXP434
#PPo
#shock=pd.read_csv("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp434/434schk/434schk_schk.csv")
#drpl=pd.read_csv("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp434/434drpl/434Gyardrp_L.csv")
#exp=pd.read_csv("/media/andres/DISPOSITIVO/mestrado/Andressss/lowpress_drpl/Mosesres/Moses434/Mosesexp434.csv")
"""
plt.figure(1)
plt.figure(figsize=(7,5))
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



xi="\u03BE"
psi="\u03C8"

plt.xlim(8.22/100,14/100)
plt.ylim(0.2,0.7)

plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

plt.plot(shock["x"]/100,shock["PPo"],"-",lw=2.5,color="darkblue",label="Weak detonation")
plt.plot(drpl["x"],drpl["PPo"],"-",lw=2.5,color="steelblue",label="Droplet growth")
plt.plot(exp["x"]/100,exp["ppo"],"--",lw=2.5,color="black",label="Moses and Stein (1978)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp434/Plots/PPo434.pdf")
"""

#Nucleation rate
"""
plt.figure(1)
plt.figure(figsize=(7,5))
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



xi="\u03BE"
psi="\u03C8"

#plt.xlim(5/100,14/100)
#plt.ylim(1,1e25)
#plt.yscale("log")


plt.ylabel("Nucleation rate J (m$^{-3}$s$^{-1}$)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

#plt.plot(shock["x"]/100,shock["PPo"],"-",lw=2.5,color="forestgreen",label="Weak detonation")
plt.plot(drpl["x"],drpl["J"],"-",lw=2.5,color="steelblue",label="Droplet growth")
#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
#plt.plot(exp["x"]/100,exp["ppo"],"--",lw=2.5,color="black",label="Moses and Stein (1978)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)


plt.legend(frameon=False,prop=font2)
plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp434/Plots/Ncl434.pdf")
"""

"""
plt.figure(1)
plt.figure(figsize=(7,5))
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



xi="\u03BE"
psi="\u03C8"

plt.xlim(10/100,16/100)
plt.ylim(0.2e-08,0.8e-08)
#plt.yscale("log")


plt.ylabel("Droplet mean radius (m)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

#plt.plot(shock["x"]/100,shock["PPo"],"-",lw=2.5,color="forestgreen",label="Weak detonation")
plt.plot(drpl["x"],drpl["rnew"],"-",lw=2.5,color="steelblue",label="Droplet growth")
#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
#plt.plot(exp["x"]/100,exp["ppo"],"--",lw=2.5,color="black",label="Moses and Stein (1978)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)
plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/ResultsMosescomp/exp434/Plots/rm434.pdf")
"""