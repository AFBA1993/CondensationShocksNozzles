#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 19:25:10 2020

@author: andres
"""

from scipy.optimize import fsolve
from Coollib import Cool_PT, Cool_Ps, Cool_hs, Cool_Psat, Cool_Tsat
from cflow_cool import fnd_isenexp, fnd_throat_MP, fnd_stagthr,fnd_throat
from Geos import new_lietteri,lietteri_expansion
from Twophase import fnd_2phase_exp, Propmix, sol2Phase
from TwophGeos import lietteri_2phase
from Cnd_shock import cshock
from solution import sol_profiles
from Rankine_line import Rayleighl, Hugoniotl, fnd_CJ, fndPs_RHeq, RayleighlS, fnd_RHS, fnd_CJS
import numpy as np
import matplotlib.pyplot as plt
from momentum import Err_mom
import pandas as pd

f="CarbonDioxide"
EoS="HEOS"


Po=57.24*100000
To=303.9
propo=Cool_PT(Po,To,f,EoS)
so=propo[4]

Pi=0.98*Po
Pf=0.9*Po
ne=35
stp=(Pf-Pi)
upT=To+30
for i in range (ne):
    prop=Cool_Ps(Pi,so,f,EoS,upT)
    print(i)
    Pi+=stp
    upT+=-0.5