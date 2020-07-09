# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 13:47:47 2020

@author: davem
"""

import matplotlib.pyplot as plt
import numpy as np

## Uncomment the following if running this script to start on its own. 
exec(open('T Cell Expansion Functions v1.2.py').read())

## This first set of following code will simulate the expansion of mouse T cells
## over the course of 7 days from a starting population of 100,000 cells, 1% of 
## which are Tregs, 16% are effector CD4, and 45.2% are CD8s.
    
[Ns, Rs, E8s, E4s, ts, params] = ExpansionSim(0.02, 0.452, 0.30)

[cell_totals, f_N, f_E8, f_E4, f_R] = FractionCalcs(Ns, Rs, E4s, E8s)

countsFig = plt.figure()
counts = plt.subplot(111)
counts.plot(ts, E8s, 'k', label = 'CD8 T Cells')
counts.plot(ts, Ns,'b', label = 'Naive CD4 T Cells')
counts.plot(ts, E4s, 'r', label = 'Effector CD4 T Cells')
counts.plot(ts, Rs, 'g', label = 'Regulatory T Cells')
counts.legend(loc = 'best')

freqsFig = plt.figure()
freqs = plt.subplot(111)
freqs.plot(ts, f_E8, 'k', label = 'CD8 T Cells')
freqs.plot(ts, f_N,'b', label = 'Naive CD4 T Cells')
freqs.plot(ts, f_E4, 'r', label = 'Effector CD4 T Cells')
freqs.plot(ts, f_R, 'g', label = 'Regulatory T Cells')
freqs.legend(loc = 'best')
#
##I'm going to save data generated from this as an np.array which just makes it easier to save
## and re-import later.
Control_Data = np.array([Ns, Rs, E4s, E8s, f_N, f_R, f_E4, f_E8, ts])

#Since most of our fractions are fractions of CD4s, I'm making an extra data storage matrix
#Where the population fractions are calculated for just the CD4 t cells. 

f_tot_CD4s = f_N + f_E4 + f_R
f_N_CD4s = f_N/f_tot_CD4s
f_R_CD4s = f_R/f_tot_CD4s
f_E4_CD4s = f_E4/f_tot_CD4s

CD4s_Control_Data = np.array([f_N_CD4s, f_R_CD4s, f_E4_CD4s, ts])

## This second set of code will simulate the the same expansion as before, but now there
## is TGF beta in the media, which I have accounted for only in its ability to enhance
## naive T cell conversion, but not in suppression of effector T cells. I also have a variable
## that accounts for Treg suppression of effector T cell expansion, but right now I have it 
## set to zero, and I might just leave it out for now.


[Ns, Rs, E8s, E4s, ts, params] = ExpansionSim(0.02, 0.452, 0.30, TGF = True)

[cell_totals, f_N, f_E8, f_E4, f_R] = FractionCalcs(Ns, Rs, E4s, E8s)

countsFig = plt.figure()
counts = plt.subplot(111)
counts.plot(ts, E8s, 'k', label = 'CD8 T Cells')
counts.plot(ts, Ns,'b', label = 'Naive CD4 T Cells')
counts.plot(ts, E4s, 'r', label = 'Effector CD4 T Cells')
counts.plot(ts, Rs, 'g', label = 'Regulatory T Cells')
counts.legend(loc = 'best')

freqsFig = plt.figure()
freqs = plt.subplot(111)
freqs.plot(ts, f_E8, 'k', label = 'CD8 T Cells')
freqs.plot(ts, f_N,'b', label = 'Naive CD4 T Cells')
freqs.plot(ts, f_E4, 'r', label = 'Effector CD4 T Cells')
freqs.plot(ts, f_R, 'g', label = 'Regulatory T Cells')
freqs.legend(loc = 'best')

TGF_Data= np.array([Ns, Rs, E4s, E8s, f_N, f_R, f_E4, f_E8, ts])

f_tot_CD4s = f_N + f_E4 + f_R
f_N_CD4s = f_N/f_tot_CD4s
f_R_CD4s = f_R/f_tot_CD4s
f_E4_CD4s = f_E4/f_tot_CD4s

CD4s_TGF_Data = np.array([f_N_CD4s, f_R_CD4s, f_E4_CD4s, ts])
# This third set of code will simulate the the same expansion as before, but now there
# is only rapamycin, which effectively slows down the expansion of effector T cells, but not
# the regulatory T cells. However, there is nothing inducing regulatory T cells in this one,
# so the population dynamics of the Tregs should be the same as the first simulation.

[Ns, Rs, E8s, E4s, ts, params] = ExpansionSim(0.02, 0.452, 0.30, Rapa = True)

[cell_totals, f_N, f_E8, f_E4, f_R] = FractionCalcs(Ns, Rs, E4s, E8s)

countsFig = plt.figure()
counts = plt.subplot(111)
counts.plot(ts, E8s, 'k', label = 'CD8 T Cells')
counts.plot(ts, Ns,'b', label = 'Naive CD4 T Cells')
counts.plot(ts, E4s, 'r', label = 'Effector CD4 T Cells')
counts.plot(ts, Rs, 'g', label = 'Regulatory T Cells')
counts.legend(loc = 'best')

freqsFig = plt.figure()
freqs = plt.subplot(111)
freqs.plot(ts, f_E8, 'k', label = 'CD8 T Cells')
freqs.plot(ts, f_N,'b', label = 'Naive CD4 T Cells')
freqs.plot(ts, f_E4, 'r', label = 'Effector CD4 T Cells')
freqs.plot(ts, f_R, 'g', label = 'Regulatory T Cells')
freqs.legend(loc = 'best')

Rapa_Data = np.array([Ns, Rs, E4s, E8s, f_N, f_R, f_E4, f_E8, ts])

f_tot_CD4s = f_N + f_E4 + f_R
f_N_CD4s = f_N/f_tot_CD4s
f_R_CD4s = f_R/f_tot_CD4s
f_E4_CD4s = f_E4/f_tot_CD4s

CD4s_Rapa_Data = np.array([f_N_CD4s, f_R_CD4s, f_E4_CD4s, ts])

# This Last set of code will simulate the expansion with both Rapa and TGF beta present
# and should take into account both suppression of effector T cell expansion due to rapa
# as well as phenotypic conversion to Tregs due to TGF.

[Ns, Rs, E8s, E4s, ts, params] = ExpansionSim(0.02, 0.452, 0.30,TGF = True, Rapa = True)

[cell_totals, f_N, f_E8, f_E4, f_R] = FractionCalcs(Ns, Rs, E4s, E8s)

countsFig = plt.figure()
counts = plt.subplot(111)
counts.plot(ts, E8s, 'k', label = 'CD8 T Cells')
counts.plot(ts, Ns,'b', label = 'Naive CD4 T Cells')
counts.plot(ts, E4s, 'r', label = 'Effector CD4 T Cells')
counts.plot(ts, Rs, 'g', label = 'Regulatory T Cells')
counts.legend(loc = 'best')

freqsFig = plt.figure()
freqs = plt.subplot(111)
freqs.plot(ts, f_E8, 'k', label = 'CD8 T Cells')
freqs.plot(ts, f_N,'b', label = 'Naive CD4 T Cells')
freqs.plot(ts, f_E4, 'r', label = 'Effector CD4 T Cells')
freqs.plot(ts, f_R, 'g', label = 'Regulatory T Cells')
freqs.legend(loc = 'best')

Combo_Data = np.array([Ns, Rs, E4s, E8s, f_N, f_R, f_E4, f_E8, ts])

f_tot_CD4s = f_N + f_E4 + f_R
f_N_CD4s = f_N/f_tot_CD4s
f_R_CD4s = f_R/f_tot_CD4s
f_E4_CD4s = f_E4/f_tot_CD4s

CD4s_Combo_Data = np.array([f_N_CD4s, f_R_CD4s, f_E4_CD4s, ts])