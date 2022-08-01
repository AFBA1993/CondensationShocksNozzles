#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 17:10:18 2020

@author: andres
"""
import matplotlib.pyplot as plt
import pandas as pd

################################
##### Mesh  ####################
################################




mesh150=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Resultscomparisonsdrpl/Mesh/gyar100moore/GyarMoore100_L.csv")
mesh200=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Resultscomparisonsdrpl/Mesh/gyar200moore/GyarMoore200_L.csv")
mesh300=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Resultscomparisonsdrpl/Mesh/gyar300moore/GyarMoore300_L.csv")
mesh400=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Resultscomparisonsdrpl/Mesh/gyar400moore/GyarMoore400_L.csv")
mesh500=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Resultscomparisonsdrpl/Mesh/gyar500moore/GyarMoore500_L.csv")
mesh1000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Resultscomparisonsdrpl/Mesh/gyar1000moore/GyarMoore1000_L.csv")
mesh2000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Resultscomparisonsdrpl/Mesh/gyar2000moore/GyarMoore2000_L.csv")
mesh3000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Resultscomparisonsdrpl/Mesh/Gyar3000moore/GyarMoore3000_L.csv")
mesh4000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Resultscomparisonsdrpl/Mesh/Gyar4000moore/GyarMoore4000_L.csv")
mesh5000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Resultscomparisonsdrpl/Mesh/Gyar5000moore/GyarMoore5000_L.csv")
mesh10000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Resultscomparisonsdrpl/Mesh/Gyar10000moore/GyarMoore10000_L.csv")




#mesh1000=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/mallascomp18c/18c1000/Gyar18cdrpl1000_L.csv")
#mesh200=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/mallascomp18c/18c200/Gyar18cdrpl200.csv_L.csv")
#mesh2500=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/mallascomp18c/18c2500/Gyar18cdrpl2500_L.csv")
#mesh5000=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/mallascomp18c/18c5000/Gyar18cdrpl5000_L.csv")
#mesh10000=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/mallascomp18c/18c10000/Gyar18cdrpl10000_L.csv")


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

plt.xlim(0.03,0.037)
plt.ylim(0.36,0.42)

plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

y="PPo"

plt.plot(mesh150["x"]/10,mesh150[y],"--", color="black",lw=1.5,label="Mesh=100")
plt.plot(mesh200["x"]/10,mesh200[y],"--", color="steelblue",lw=1.5,label="Mesh=200")
plt.plot(mesh300["x"]/10,mesh300[y],"--", color="darkred",lw=1.5,label="Mesh=300")
plt.plot(mesh400["x"]/10,mesh400["PPo"],"-",lw=1.0,color="darkgreen",label="Mesh=400")
#plt.plot(mesh500["x"]/10,mesh500[y],"-",lw=1.0,color="forestgreen",label="Mesh=500")
plt.plot(mesh1000["x"]/10,mesh1000[y],"-",color="dodgerblue",lw=1.2,label="Mesh=1000")
#plt.plot(mesh2000["x"]/10,mesh2000[y],"--",color="green",lw=1.2,label="Mesh=2000")
#plt.plot(mesh3000["x"],mesh3000["PPo"],"-",color="darkgreen",lw=1.2,label="Mesh=3000")
plt.plot(mesh4000["x"]/10,mesh4000[y],"-",color="darkcyan",lw=1.2,label="Mesh=4000")
plt.plot(mesh5000["x"]/10,mesh5000[y],"-",color="darkblue",lw=1.0,label="Mesh=5000")
plt.plot(mesh10000["x"]/10,mesh10000[y],"-",color="forestgreen",lw=1.0,label="Mesh=10000")

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.savefig("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/Resultscomparisonsdrpl/Plotsresults/meshmoore.pdf")


















"""
young62alpha58=pd.read_csv("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/Resultscomparisons/ResultsMoore/Young62alpha58/Young62alpha58_L.csv")
young62=pd.read_csv("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/Resultscomparisons/ResultsMoore/Young62alpha4/Young62alpha4_L.csv")
Gyarmathy60=pd.read_csv("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/Resultscomparisons/ResultsMoore/Gyarmathy60/Gyarmathy60_L.csv")
Yang=pd.read_csv("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/Resultscomparisons/ResultsMoore/Yang615/Yang615_L.csv")


exp=pd.read_csv("/media/andres/DISPOSITIVO/mestrado/Andressss/lowpress_drpl/Mooreres/MEXP.csv")
expr=pd.read_csv("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/Mooreradius.csv")
"""




"s"
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

plt.xlim(0.2,0.6)
plt.ylim(0.2,0.6)

plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

plt.plot(young62alpha58["x"],young62alpha58["PPo"],"--",lw=2.5,color="forestgreen",label="Young "+xi+"=0.62 "+psi+"=5.8")
plt.plot(Yang["x"],Yang["PPo"],"-",lw=2.5,color="dodgerblue",label="Yang & Shen "+xi+"=0.615")
plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
plt.plot(exp["x"]/10,exp["ppo"],"o",color="black",label="Moore et al. (1973)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
#plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/Resultscomparisons/Plotsresults/PPocomp.pdf")
"""
####################################################################################################################

#Nucleation rate

"""
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

plt.xlim(0.15,0.6)
plt.ylim(1,1e23)
plt.yscale("log")

plt.ylabel("Nucleation rate J (m$^{-3}$s$^{-1}$)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

plt.plot(young62alpha58["x"],young62alpha58["J"],"--",lw=2.5,color="forestgreen",label="Young "+xi+"=0.62 "+psi+"=5.8")
plt.plot(Yang["x"],Yang["J"],"-",lw=2.5,color="dodgerblue",label="Yang & Shen "+xi+"=0.615")
plt.plot(Gyarmathy60["x"],Gyarmathy60["J"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/Resultscomparisons/Plotsresults/Jcomplog.pdf")
"""
###########################################################################################################
#Number of droplets

"""
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

plt.xlim(0.1,0.7)
plt.ylim(1,2e20)
plt.yscale("log")

plt.ylabel("Number of droplets N (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

plt.plot(young62alpha58["x"],young62alpha58["N"],"--",lw=2.5,color="forestgreen",label="Young "+xi+"=0.62 "+psi+"=5.8")
plt.plot(Yang["x"],Yang["N"],"-",lw=2.5,color="dodgerblue",label="Yang & Shen "+xi+"=0.615")
plt.plot(Gyarmathy60["x"],Gyarmathy60["N"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/Resultscomparisons/Plotsresults/NcompS.pdf")
"""
####################################################################################################################
#Wetness mass fraction

"""
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

#plt.xlim(0.1,0.7)
#plt.ylim(1,2e20)
#plt.yscale("log")

plt.ylabel("Wetness fraction y (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

plt.plot(young62alpha58["x"],young62alpha58["y"],"--",lw=2.5,color="forestgreen",label="Young "+xi+"=0.62 "+psi+"=5.8")
plt.plot(Yang["x"],Yang["y"],"-",lw=2.5,color="dodgerblue",label="Yang & Shen "+xi+"=0.615")
plt.plot(Gyarmathy60["x"],Gyarmathy60["y"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/Resultscomparisons/Plotsresults/ycomp.pdf")

"""

#########################################################################################################
#Supersaturation
"""
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

#plt.xlim(0.1,0.7)
#plt.ylim(1,2e20)
#plt.yscale("log")

plt.ylabel("Supersaturation ratio  S (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

plt.plot(young62alpha58["x"],young62alpha58["S"],"--",lw=2.5,color="forestgreen",label="Young "+xi+"=0.62 "+psi+"=5.8")
plt.plot(Yang["x"],Yang["S"],"-",lw=2.5,color="dodgerblue",label="Yang & Shen "+xi+"=0.615")
plt.plot(Gyarmathy60["x"],Gyarmathy60["S"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/Resultscomparisons/Plotsresults/Scomp.pdf")
"""
#####################################
#Subcooling
"""
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

#plt.xlim(0.1,0.7)
#plt.ylim(1,2e20)
#plt.yscale("log")

plt.ylabel("Subcooling $\Delta$T (K)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

plt.plot(young62alpha58["x"],young62alpha58["DelT"],"--",lw=2.5,color="forestgreen",label="Young "+xi+"=0.62 "+psi+"=5.8")
plt.plot(Yang["x"],Yang["DelT"],"-",lw=2.5,color="dodgerblue",label="Yang & Shen "+xi+"=0.615")
plt.plot(Gyarmathy60["x"],Gyarmathy60["DelT"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/Resultscomparisons/Plotsresults/Subcooling.pdf")
"""
###############################################################
#Radius

"""
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

plt.xlim(0.25,0.7)
plt.ylim(0.4E-8,7e-8)
#plt.yscale("log")

plt.ylabel("Droplet mean radius r (m)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

plt.plot(young62alpha58["x"],young62alpha58["rnew"],"--",lw=2.5,color="forestgreen",label="Young "+xi+"=0.62 "+psi+"=5.8")
plt.plot(Yang["x"],Yang["rnew"],"-",lw=2.5,color="dodgerblue",label="Yang & Shen "+xi+"=0.615")
plt.plot(Gyarmathy60["x"],Gyarmathy60["rnew"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
plt.plot(expr["x"]/10,expr["rm"],"o",color="black",label="Moore et al. (1973)")


plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/Resultscomparisons/Plotsresults/meanradius.pdf")
"""
######################################


"""
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

#plt.xlim(0.25,0.7)
#plt.ylim(0.4E-8,7e-8)
#plt.yscale("log")

plt.ylabel("Specific entropy s (J/kg K )",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

plt.plot(young62alpha58["x"],young62alpha58["s"],"--",lw=2.5,color="forestgreen",label="Young "+xi+"=0.62 "+psi+"=5.8")
plt.plot(Yang["x"],Yang["s"],"-",lw=2.5,color="dodgerblue",label="Yang & Shen "+xi+"=0.615")
plt.plot(Gyarmathy60["x"],Gyarmathy60["s"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
#plt.plot(expr["x"]/10,expr["rm"],"o",color="black",label="Moore et al. (1973)")


plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.savefig("/media/andres/DISPOSITIVO/mestrado/Andressss/Final/Resultscomparisons/Plotsresults/entropy.pdf")
"""


