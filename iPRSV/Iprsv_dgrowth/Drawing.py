#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 17:10:18 2020

@author: andres
"""
import matplotlib.pyplot as plt
import pandas as pd


############################################################
######################### IPRSV-MOORE ####################
############################################################




mesh150=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moore_iprsv/Mesh_iprsv_Moore/iprsvmoore150/dprliprsvM150_L.csv")
mesh200=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moore_iprsv/Mesh_iprsv_Moore/iprsvmooore200/dprliprsvM200_L.csv")
mesh300=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moore_iprsv/Mesh_iprsv_Moore/iprsvmoore300/dprliprsvM300_L.csv")
mesh400=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moore_iprsv/Mesh_iprsv_Moore/iprsvmoore400/dprliprsvM400_L.csv")
mesh500=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moore_iprsv/Mesh_iprsv_Moore/iprsvmoore500/dprliprsvM500_L.csv")
mesh1000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moore_iprsv/Mesh_iprsv_Moore/iprsvmooore1000/dprliprsvM1000_L.csv")
mesh2000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moore_iprsv/Mesh_iprsv_Moore/iprsvmoore2000/dprliprsvM2000_L.csv")
mesh3000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moore_iprsv/Mesh_iprsv_Moore/iprsvmoore3000/iprsvmooore4000dprliprsvM3000_L.csv")
mesh4000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moore_iprsv/Mesh_iprsv_Moore/iprsvmooore4000/iprsvmooore4000dprliprsvM4000_L.csv")
mesh5000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moore_iprsv/Mesh_iprsv_Moore/iprsvmooore5000/iprsvmooore5000dprliprsvM5000_L.csv")

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

plt.xlim(2.7/10,3.4/10)
plt.ylim(0.36,0.44)

plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)


plt.plot(mesh150["x"]/10,mesh150["PPo"],"--", color="blue",lw=1.5,label="Mesh=150")
plt.plot(mesh200["x"]/10,mesh200["PPo"],"--", color="steelblue",lw=1.5,label="Mesh=200")
plt.plot(mesh300["x"]/10,mesh300["PPo"],"--", color="darkred",lw=1.5,label="Mesh=1000")
#plt.plot(mesh400["x"],mesh400["PPo"],"-",lw=1.0,color="black",label="Mesh=2500")
plt.plot(mesh500["x"]/10,mesh500["PPo"],"-",lw=1.0,color="forestgreen",label="Mesh=2500")
plt.plot(mesh1000["x"]/10,mesh1000["PPo"],"--",color="dodgerblue",lw=1.2,label="Mesh=1000")
plt.plot(mesh2000["x"]/10,mesh2000["PPo"],"--",color="green",lw=1.2,label="Mesh=2000")
#plt.plot(mesh3000["x"],mesh3000["PPo"],"-",color="darkgreen",lw=1.2,label="Mesh=3000")
plt.plot(mesh4000["x"]/10,mesh4000["PPo"],"-",color="darkcyan",lw=1.2,label="Mesh=4000")
plt.plot(mesh5000["x"]/10,mesh5000["PPo"],"-",color="darkblue",lw=1.0,label="Mesh=5000")

"""
plt.plot(mesh5000["x"],mesh5000["PPo"],"--",color="dodgerblue",lw=1.2,label="Mesh=5000")
plt.plot(mesh10000["x"],mesh10000["PPo"],"-",color="green",lw=1.0,label="Mesh=10000")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
"""
"""
plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

#plt.savefig("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/IPRSV_RESULTS/Moore_iprsv/Mesh_iprsv_Moore/plotsiprsv/iprsvmoore.pdf")
"""






"""
schk=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/IPRSV_RESULTS/Moore_iprsv/Moore_shck/iprsvmoore_schk.csv")
drpl=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/IPRSV_RESULTS/Moore_iprsv/Mesh_iprsv_Moore/iprsvmooore4000/iprsvmooore4000dprliprsvM4000_L.csv")
exp=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/lowpress_drpl/Mooreres/MEXP.csv")
#expr=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Higpress_drpl/res18c/expr_18c.csv")
"""

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

plt.xlim(0.2,0.6)
plt.ylim(0.2,0.7)

plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)

plt.plot(schk["x"]/10,schk["PPo"],"-",lw=2.5,color="dodgerblue",label="Weak detonation")
plt.plot(drpl["x"]/10,drpl["PPo"],"-",lw=2.5,color="darkcyan",label="Droplet growth")#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
plt.plot(exp["x"]/10,exp["ppo"],"o",lw=2.5,color="black",label="Moore et al. (1973)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)

plt.savefig("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/IPRSV_RESULTS/Moore_iprsv/Mesh_iprsv_Moore/plotsiprsv/iprsvmooreppo.pdf")

#radious
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

plt.xlim(2.77/10,0.7)
plt.ylim(1e-09,6.5e-09)

plt.ylabel("Droplet mean radius (m)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)


plt.plot(drpl["x"]/10,drpl["rm"],"-",lw=2.5,color="darkcyan",label="Droplet growth")#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
#plt.plot(expr["x"],expr["r"],"o",lw=2.5,color="black",label="Gyarmathy (2005)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.savefig("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/IPRSV_RESULTS/Moore_iprsv/Mesh_iprsv_Moore/plotsiprsv/iprsvmoorerm.pdf")
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

plt.xlim(0.225,0.350)
#plt.ylim(0.5e-08,7e-08)

plt.ylabel("Nucleation rate J (m$^{-3}$s$^{-1}$)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (m)",fontsize=lsize,family=family)


plt.plot(drpl["x"]/10,drpl["J"],"-",lw=2.5,color="darkcyan",label="Droplet growth")#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
#plt.plot(exp["x"],exp["J"],"o",lw=2.5,color="black",label="Gyarmathy (2005)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
#plt.savefig("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/IPRSV_RESULTS/Moore_iprsv/Mesh_iprsv_Moore/plotsiprsv/nrate.pdf")

"""

############################################################
######################### IPRSV-moses ####################
############################################################

mesh150=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moses_iprsv/moses_iprsv100/mosesprliprsv100_L.csv")
mesh200=pd.read_csv("//media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moses_iprsv/moses_iprsv200/mosesprliprsv200_L.csv")
mesh300=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moses_iprsv/moses_iprsv300/mosesprliprsv300_L.csv")
mesh400=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moses_iprsv/moses_iprsv400/mosesprliprsv400_L.csv")
mesh500=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moses_iprsv/moses_iprsv500/mosesprliprsv500_L.csv")
mesh1000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moses_iprsv/moses_iprsv1000/mosesprliprsv1000_L.csv")
mesh2000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moses_iprsv/moses_iprsv2000/mosesprliprsv2000_L.csv")
mesh3000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moses_iprsv/moses_iprsv3000/mosesprliprsv3000_L.csv")
mesh4000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moses_iprsv/mosesiprsv4000/mosesprliprsv400_L.csv")
mesh5000=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moses_iprsv/mosesiprsv5000/mosesprliprsv5000_L.csv")

#mesh1000=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/mallascomp18c/18c1000/Gyar18cdrpl1000_L.csv")
#mesh200=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/mallascomp18c/18c200/Gyar18cdrpl200.csv_L.csv")
#mesh2500=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/mallascomp18c/18c2500/Gyar18cdrpl2500_L.csv")
#mesh5000=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/mallascomp18c/18c5000/Gyar18cdrpl5000_L.csv")
#mesh10000=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Final/ResultsHIGH/Result18c/mallascomp18c/18c10000/Gyar18cdrpl10000_L.csv")

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



xi="\u03BE"
psi="\u03C8"

#plt.xlim(10.0,11.5)
#plt.ylim(0.35,0.43)

plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (cm)",fontsize=lsize,family=family)
y="J"

plt.plot(mesh150["x"],mesh150[y],"--", color="blue",lw=1.5,label="Mesh=100")
plt.plot(mesh200["x"],mesh200[y],"--", color="steelblue",lw=1.5,label="Mesh=200")
#plt.plot(mesh300["x"]/10,mesh300["PPo"],"--", color="darkred",lw=1.5,label="Mesh=1000")
#plt.plot(mesh400["x"],mesh400["PPo"],"-",lw=1.0,color="black",label="Mesh=2500")
plt.plot(mesh500["x"],mesh500[y],"-",lw=1.0,color="forestgreen",label="Mesh=500")
#plt.plot(mesh1000["x"],mesh1000["PPo"],"--",color="darkred",lw=1.2,label="Mesh=1000")
plt.plot(mesh2000["x"],mesh2000[y],"--",color="dodgerblue",lw=1.2,label="Mesh=2000")
plt.plot(mesh3000["x"],mesh3000[y],"-",color="darkred",lw=1.2,label="Mesh=3000")
plt.plot(mesh4000["x"],mesh4000[y],"-",color="darkcyan",lw=1.2,label="Mesh=4000")
plt.plot(mesh5000["x"],mesh5000[y],"-",color="darkblue",lw=1.0,label="Mesh=5000")


plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
#plt.savefig("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moses_iprsv/mosesplots/meshmoses417.pdf")

"""





















schk=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moses_iprsv/moseschk/iprsvmoses417_schk.csv")
drpl=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moses_iprsv/mosesdprl/mosesdprliprsv4000_L.csv")
exp=pd.read_csv("/media/andres777/DISPOSITIVO/Mestrado/Andressss/lowpress_drpl/Mosesres/Moses417/Mosesexp417.csv")
#expr=pd.read_csv("/home/andres/Área de Trabalho/BACKUP/Mestrado/Andressss/Higpress_drpl/res18c/expr_18c.csv")


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

plt.xlim(8.22,14)
plt.ylim(0.1,0.7)

plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (cm)",fontsize=lsize,family=family)

plt.plot(schk["x"],schk["PPo"],"-",lw=2.5,color="slategrey",label="Weak detonation")
plt.plot(drpl["x"],drpl["PPo"],"-",lw=2.5,color="darkcyan",label="Droplet growth")#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
plt.plot(exp["x"],exp["ppo"],"--",lw=2.5,color="black",label="Moses and Stein (1978)")
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)


plt.savefig("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moses_iprsv/mosesplots/ppo417.pdf")
"""




#radius

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

plt.xlim(10.3,16)
plt.ylim(1e-09,0.7e-08)

plt.ylabel("Droplet mean radius (m)",fontsize=lsize,family=family)
#plt.ylabel("Pressure ratio P/Po (-)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (cm)",fontsize=lsize,family=family)

#plt.plot(schk["x"],schk["PPo"],"-",lw=2.5,color="slategrey",label="Weak detonation")
plt.plot(drpl["x"],drpl["rm"],"-",lw=2.5,color="darkcyan",label="Droplet growth")#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
#plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.savefig("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moses_iprsv/mosesplots/RM.pdf")
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

plt.xlim(9.0,12)
#plt.ylim(1e-09,0.7e-08)

#plt.ylabel("Droplet mean radius (m)",fontsize=lsize,family=family)
plt.ylabel("Nucleation rate J (m$^{-3}$s$^{-1}$)",fontsize=lsize,family=family)
plt.xlabel("Nozzle axis (cm)",fontsize=lsize,family=family)

#plt.plot(schk["x"],schk["PPo"],"-",lw=2.5,color="slategrey",label="Weak detonation")
plt.plot(drpl["x"],drpl["J"],"-",lw=2.5,color="darkcyan",label="Droplet growth")#plt.plot(Gyarmathy60["x"],Gyarmathy60["PPo"],"--",lw=2.5,color="navy",label="Gyarmathy "+xi+"=0.6")
#plt.grid(color="gainsboro", linestyle='--', linewidth=.75)
plt.grid(color="gainsboro", linestyle='--', linewidth=.75)

plt.xticks(fontsize=axisize,family=family)
plt.yticks(fontsize=axisize,family=family)

plt.legend(frameon=False,prop=font2)
plt.savefig("/media/andres777/DISPOSITIVO/Mestrado/Andressss/Final/IPRSV_RESULTS/Moses_iprsv/mosesplots/nrate.pdf")
"""








