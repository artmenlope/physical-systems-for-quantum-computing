# -*- coding: utf-8 -*-
"""
@author: Arturo Mena López

Main program for calling transmon_and_cpb.py as a
module (must be in the same directory).
"""

from transmon_and_cpb import *
import numpy as np
import matplotlib.pyplot as plt


#%%
# =============================================================================
# COOPER PAIR BOX
# =============================================================================

# Parameters.
N  = 100 # Number of spatial points (for the ham() function in transmon_and_cpb.py).
Ej = 1 # Josephson energy in GHz (hbar=1).
Ec = 10 # Charging energy in GHz (hbar=1).
kval = 3 # Number of plotted energy states.
phimax = np.pi # Maximum value of the phase variable in the plots. The minimum will be -phimax.
ng_max = 1 # Maximum value of the charge offset in the plots. The minimum will be -ng_max.
num_ngs = N # Number of points in the ng-space.

enerlist, EmE01_list = plot_Em_ng(ng_max=ng_max, kval=kval, Ej=Ej, Ec=Ec, num_ngs=500, N=N, phimax=phimax)


#%%
# =============================================================================
# TRANSMON
# =============================================================================

# Parameters.
N  = 500 # Number of spatial points (for the ham() function in transmon_and_cpb.py).
Ej = 50 # Josephson energy in GHz (hbar=1).
Ec = 1 # Charging energy in GHz (hbar=1).
ng = 1/2 # Charge offset (for the transmon).
phimax = np.pi # Maximum value of the phase variable in the plots. The minimum will be -phimax.
ng_max = 2 # Maximum value of the charge offset in the plots. The minimum will be -ng_max.
num_ngs = N # Number of points in the ng-space.


kval = 6 # Number of plotted energy states.
mat = ham(phimax, N, potfun=josephson, param=Ej, ng=ng, Ej=Ej, Ec=Ec)
plotstates(mat,kval,phimax,N,potfun=josephson,param=Ej, Ej=Ej, Ec=Ec)


kval = 3 # Number of plotted energy states.
Ej_list = [1, 5, 10, 50] # List of Josephson energies in GHz (hbar=1).
Ec_list = [1, 1, 1, 1] # List of charging energies in GHz (hbar=1).
plot_Em_ng_subplots(ng_max=ng_max, kval=kval, Ej_list=Ej_list, Ec_list=Ec_list, num_ngs=num_ngs, N=N, phimax=phimax)
