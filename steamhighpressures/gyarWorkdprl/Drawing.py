#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 17:10:18 2020

@author: andres
"""
import matplotlib.pyplot as plt
import pandas as pd

########################################################################################################
############### Mesh comparison 18c ##################################
#############################################################

mesh100=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/mallascomp18c/18c100/Gyar18cdrpl100_L.csv")
mesh1000=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/mallascomp18c/18c1000/Gyar18cdrpl1000_L.csv")
mesh200=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/mallascomp18c/18c200/Gyar18cdrpl200.csv_L.csv")
mesh2500=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/mallascomp18c/18c2500/Gyar18cdrpl2500_L.csv")
mesh5000=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/mallascomp18c/18c5000/Gyar18cdrpl5000_L.csv")
mesh10000=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/mallascomp18c/18c10000/Gyar18cdrpl10000_L.csv")


#Pressure profile


plt.figure(1)
plt.figure(figsize=(7.5,7.5))
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

plt.xlim(4.3,7)
plt.ylim(0.4,0.44)

plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (mm)",fontsize=lsize,family=family)


plt.plot(mesh100["x"],mesh100["PPo"],"--", color="blue",lw=1.5,label="Mesh=100")
plt.plot(mesh200["x"],mesh200["PPo"],"--", color="steelblue",lw=1.5,label="Mesh=200")
plt.plot(mesh1000["x"],mesh1000["PPo"],"--", color="darkred",lw=1.5,label="Mesh=1000")
plt.plot(mesh2500["x"],mesh2500["PPo"],"-",lw=1.0,color="black",label="Mesh=2500")

plt.plot(mesh5000["x"],mesh5000["PPo"],"--",color="dodgerblue",lw=1.2,label="Mesh=5000")
plt.plot(mesh10000["x"],mesh10000["PPo"],"-",color="green",lw=1.0,label="Mesh=10000")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)

plt.savefig("/media/andres/DISPOSITIVO/Mestrado/aquivamosenlaescrita/Tesis maestria3/fig/Results/steamhigh/18c/mesh.pdf")




############################################################
######################### Gyarmathy 18c ####################
############################################################


schk=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/18cshk/Gyarmth18cshck6000.csv_schk.csv")
drpl=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/18cdprl/Gyar18cdrpl105.csv_L.csv")
exp=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Higpress_drpl/res18c/exp18c.csv")
expr=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Higpress_drpl/res18c/expr_18c.csv")


##Ppo
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

plt.xlim(0,17)
plt.ylim(0.2,0.7)

plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (mm)",fontsize=lsize,family=family)

plt.plot(schk["x"],schk["PPo"],"-",lw=2.5,color="dodgerblue",label="Weak detonation")
plt.plot(drpl["x"],drpl["PPo"],"-",lw=2.5,color="darkcyan",label="Droplet growth")#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
plt.plot(exp["x"],exp["ppo"],"o",lw=2.5,color="black",label="Gyarmathy (2005)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)

plt.savefig("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/Plots18c/PPocomp18c.pdf")
"""
#radious

"""
plt.figure(2)
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

plt.xlim(0,14)
plt.ylim(1,4E+24)

plt.ylabel("Nucleation rate J (m$^{-3}$s$^{-1}$)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (mm)",fontsize=lsize,family=family)

#plt.plot(schk["x"],schk["y"],"-",lw=2.5,color="dodgerblue",label="Weak detonation")
plt.plot(drpl["x"],drpl["J"],"-",lw=2.5,color="darkcyan",label="Droplet growth")#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
#plt.plot(exp["x"],exp["ppo"],"o",lw=2.5,color="black",label="Gyarmathy (2005)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)

plt.savefig("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/Plots18c/Nrate.pdf")
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

plt.xlim(0,30)
plt.ylim(1e-08,7e-08)

plt.ylabel("Droplet mean radius (m)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (mm)",fontsize=lsize,family=family)


plt.plot(drpl["x"],drpl["rm"],"-",lw=2.5,color="darkcyan",label="Droplet growth")#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
plt.plot(expr["x"],expr["r"],"o",lw=2.5,color="black",label="Gyarmathy (2005)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.savefig("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/Plots18c/rm18c.pdf")
"""

"""

plt.figure(2)
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

plt.xlim(0,10)
#plt.ylim(0.5e-08,7e-08)

plt.ylabel("Nucleation rate J (m$^{-3}$s$^{-1}$)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (mm)",fontsize=lsize,family=family)


plt.plot(drpl["x"],drpl["J"],"-",lw=2.5,color="darkcyan",label="Droplet growth")#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
#plt.plot(exp["x"],exp["J"],"o",lw=2.5,color="black",label="Gyarmathy (2005)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.savefig("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/Plots18c/Nr18c.pdf")
"""


############################################################
######################### Gyarmathy 18c ####################
############################################################

"""
schk=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18b/Gyar18bschk/Gyarmth18b_schk.csv")
drpl=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18b/Gyar18bdrpl/Gyar18bdrpl3_L.csv")
exp=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Higpress_drpl/res18b/exp18b.csv")
expr=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Higpress_drpl/res18b/expr_18b.csv")
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

plt.xlim(0,19)
plt.ylim(0.15,0.7)

plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (mm)",fontsize=lsize,family=family)

plt.plot(schk["x"],schk["PPo"],"-",lw=2.5,color="darkslategrey",label="Weak detonation")
plt.plot(drpl["x"],drpl["PPo"],"-",lw=2.5,color="dodgerblue",label="Droplet growth")#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
plt.plot(exp["x"],exp["ppo"],"o",lw=2.5,color="black",label="Gyarmathy (2005)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)

plt.savefig("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18b/Plots/18bPPo.pdf")
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

plt.ylim(1e-08,7e-08)
plt.xlim(5,30)

plt.ylabel("Droplet mean radius (m)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (mm)",fontsize=lsize,family=family)


plt.plot(drpl["x"],drpl["rnew"],"-",lw=2.5,color="dodgerblue",label="Droplet growth")
plt.plot(expr["x"],expr["r"],"o",lw=2.5,color="black",label="Gyarmathy (2005)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.legend(frameon=False,prop=font2)
plt.savefig("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18b/Plots18b/18brm.pdf")
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

plt.ylim(1,3.5E+24)
plt.xlim(5,15)

plt.ylabel("Nucleation rate J (m$^{-3}$s$^{-1}$)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (mm)",fontsize=lsize,family=family)


plt.plot(drpl["x"],drpl["J"],"-",lw=2.5,color="dodgerblue",label="Droplet growth")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.legend(frameon=False,prop=font2)
plt.savefig("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18b/Plots18b/18bNrate.pdf")
"""