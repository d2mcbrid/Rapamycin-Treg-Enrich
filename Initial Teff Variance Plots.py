# -*- coding: utf-8 -*-
"""
Created on Thu Apr 02 14:37:18 2020

@author: davem
"""
import matplotlib.pyplot as plt
import numpy as np
#%%
#Make sure that all fonts are set to Arial
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
#%%
#Load the data to be graphed from numpy arrays saved from previous experiments
#Rapa data loading
CD4s_Rapa_Data_m_30pctEffector = np.load('C:/Users/davem/OneDrive/Desktop/UCSD Expts/Computational Studies/CRC and TGF-Beta Model/Data Storage/CD4s_Rapa_Data_m_30pctEffector.npy')
CD4s_Rapa_Data_m_20pctEffector = np.load('C:/Users/davem/OneDrive/Desktop/UCSD Expts/Computational Studies/CRC and TGF-Beta Model/Data Storage/CD4s_Rapa_Data_m_20pctEffector.npy')
CD4s_Rapa_Data_m_10pctEffector = np.load('C:/Users/davem/OneDrive/Desktop/UCSD Expts/Computational Studies/CRC and TGF-Beta Model/Data Storage/CD4s_Rapa_Data_m_10pctEffector.npy')
#TGF loading
CD4s_TGF_Data_m_30pctEffector = np.load('C:/Users/davem/OneDrive/Desktop/UCSD Expts/Computational Studies/CRC and TGF-Beta Model/Data Storage/CD4s_TGF_Data_m_30pctEffector.npy')
CD4s_TGF_Data_m_20pctEffector = np.load('C:/Users/davem/OneDrive/Desktop/UCSD Expts/Computational Studies/CRC and TGF-Beta Model/Data Storage/CD4s_TGF_Data_m_20pctEffector.npy')
CD4s_TGF_Data_m_10pctEffector = np.load('C:/Users/davem/OneDrive/Desktop/UCSD Expts/Computational Studies/CRC and TGF-Beta Model/Data Storage/CD4s_TGF_Data_m_10pctEffector.npy')
#Combo loading
CD4s_Combo_Data_m_30pctEffector = np.load('C:/Users/davem/OneDrive/Desktop/UCSD Expts/Computational Studies/CRC and TGF-Beta Model/Data Storage/CD4s_Combo_Data_m_30pctEffector.npy')
CD4s_Combo_Data_m_20pctEffector = np.load('C:/Users/davem/OneDrive/Desktop/UCSD Expts/Computational Studies/CRC and TGF-Beta Model/Data Storage/CD4s_Combo_Data_m_20pctEffector.npy')
CD4s_Combo_Data_m_10pctEffector = np.load('C:/Users/davem/OneDrive/Desktop/UCSD Expts/Computational Studies/CRC and TGF-Beta Model/Data Storage/CD4s_Combo_Data_m_10pctEffector.npy')
#Control loading
CD4s_Control_Data_m_30pctEffector = np.load('C:/Users/davem/OneDrive/Desktop/UCSD Expts/Computational Studies/CRC and TGF-Beta Model/Data Storage/CD4s_Control_Data_m_30pctEffector.npy')
CD4s_Control_Data_m_20pctEffector = np.load('C:/Users/davem/OneDrive/Desktop/UCSD Expts/Computational Studies/CRC and TGF-Beta Model/Data Storage/CD4s_Control_Data_m_20pctEffector.npy')
CD4s_Control_Data_m_10pctEffector = np.load('C:/Users/davem/OneDrive/Desktop/UCSD Expts/Computational Studies/CRC and TGF-Beta Model/Data Storage/CD4s_Control_Data_m_10pctEffector.npy')


ts = CD4s_TGF_Data_m_20pctEffector[3,:]
#%% Code to plot variance in results based on initial conditions for Rapa
# Make our figure to call
RapaFig = plt.figure()
# Make our subplot/axes. I still need to figure out what the differences are between a subplot and an axis set
Rapa = plt.subplot(111)
# Plot the regulatory T cell data
Rapa.plot(ts, CD4s_Rapa_Data_m_30pctEffector[1,:], 'g--', label = '30% Teff')
Rapa.plot(ts, CD4s_Rapa_Data_m_20pctEffector[1,:], 'g', label = '20% Teff')
Rapa.plot(ts, CD4s_Rapa_Data_m_10pctEffector[1,:], 'g:', label = '10% Teff')
# This bit of code does the filling between the upper and lower limits to make it look pretty. 
Rapa.fill_between(ts, CD4s_Rapa_Data_m_30pctEffector[1,:],CD4s_Rapa_Data_m_10pctEffector[1,:],color = 'green', alpha = 0.2)
#Same shit but for the effector T cell data
Rapa.plot(ts, CD4s_Rapa_Data_m_30pctEffector[2,:], 'r--', label = '30% Teff')
Rapa.plot(ts, CD4s_Rapa_Data_m_20pctEffector[2,:], 'r', label = '20% Teff')
Rapa.plot(ts, CD4s_Rapa_Data_m_10pctEffector[2,:], 'r:', label = '10% Teff')
Rapa.fill_between(ts, CD4s_Rapa_Data_m_30pctEffector[2,:],CD4s_Rapa_Data_m_10pctEffector[2,:],color = 'red', alpha = 0.2)
# Make a legend
# Rapa.legend(loc = 'best', fontsize = 16)
#Add axis titles
Rapa.set_xlabel('Time (Hours)', size = 22)
Rapa.set_ylabel('Fraction of CD4 T cells', size = 22)

#Set axis sizes and where the tick marks are
Rapa.set_xticks(ticks = [0,50,100,150])
Rapa.set_xticklabels(labels = [0,50,100,150], fontsize = 18)
Rapa.set_yticks(ticks = [0,0.2,0.4,0.6,0.8,1.0])
Rapa.set_yticklabels(labels = [0,0.2,0.4,0.6,0.8,1.0], fontsize = 18)
#Set how far the axis goes
Rapa.set_ylim(0,1.0)
Rapa.set_xlim(0,150)
RapaFig.tight_layout()
RapaFig.savefig('Rapa_Teff_variance', dpi = 1200)
#%% Code to plot variance in results based on initial conditions for TGF
TGFFig = plt.figure()
TGF = plt.subplot(111)
TGF.plot(ts, CD4s_TGF_Data_m_30pctEffector[1,:], 'g--', label = '30% Teff')
TGF.plot(ts, CD4s_TGF_Data_m_20pctEffector[1,:], 'g', label = '20% Teff')
TGF.plot(ts, CD4s_TGF_Data_m_10pctEffector[1,:], 'g:', label = '10% Teff')
TGF.fill_between(ts, CD4s_TGF_Data_m_30pctEffector[1,:],CD4s_TGF_Data_m_10pctEffector[1,:],color = 'green', alpha = 0.2)
TGF.plot(ts, CD4s_TGF_Data_m_30pctEffector[2,:], 'r--', label = '30% Teff')
TGF.plot(ts, CD4s_TGF_Data_m_20pctEffector[2,:], 'r', label = '20% Teff')
TGF.plot(ts, CD4s_TGF_Data_m_10pctEffector[2,:], 'r:', label = '10% Teff')
TGF.fill_between(ts, CD4s_TGF_Data_m_30pctEffector[2,:],CD4s_TGF_Data_m_10pctEffector[2,:],color = 'red', alpha = 0.2)
TGF.set_xlabel('Time (Hours)', size = 22)
TGF.set_ylabel('Fraction of CD4 T cells', size = 22)
TGF.set_xticks(ticks = [0,50,100,150])
TGF.set_xticklabels(labels = [0,50,100,150], fontsize = 18)
TGF.set_yticks(ticks = [0,0.2,0.4,0.6,0.8,1.0])
TGF.set_yticklabels(labels = [0,0.2,0.4,0.6,0.8,1.0], fontsize = 18)
TGF.set_ylim(0,1.0)
TGF.set_xlim(0,150)
TGFFig.tight_layout()
TGFFig.savefig('TGF_Teff_variance',dpi = 1200)
#%% Code to plot variance in results based on initial conditions
ComboFig = plt.figure()
Combo = plt.subplot(111)
Combo.plot(ts, CD4s_Combo_Data_m_30pctEffector[1,:], 'g--', label = '30% Teff')
Combo.plot(ts, CD4s_Combo_Data_m_20pctEffector[1,:], 'g', label = '20% Teff')
Combo.plot(ts, CD4s_Combo_Data_m_10pctEffector[1,:], 'g:', label = '10% Teff')
Combo.fill_between(ts, CD4s_Combo_Data_m_30pctEffector[1,:],CD4s_Combo_Data_m_10pctEffector[1,:],color = 'green', alpha = 0.2)
Combo.plot(ts, CD4s_Combo_Data_m_30pctEffector[2,:], 'r--', label = '30% Teff')
Combo.plot(ts, CD4s_Combo_Data_m_20pctEffector[2,:], 'r', label = '20% Teff')
Combo.plot(ts, CD4s_Combo_Data_m_10pctEffector[2,:], 'r:', label = '10% Teff')
Combo.fill_between(ts, CD4s_Combo_Data_m_30pctEffector[2,:],CD4s_Combo_Data_m_10pctEffector[2,:],color = 'red', alpha = 0.2)
Combo.set_xlabel('Time (Hours)', size = 22)
Combo.set_ylabel('Fraction of CD4 T cells', size = 22)
Combo.set_xticks(ticks = [0,50,100,150])
Combo.set_xticklabels(labels = [0,50,100,150], fontsize = 18)
Combo.set_yticks(ticks = [0,0.2,0.4,0.6,0.8,1.0])
Combo.set_yticklabels(labels = [0,0.2,0.4,0.6,0.8,1.0], fontsize = 18)
Combo.set_ylim(0,1.0)
Combo.set_xlim(0,150)
ComboFig.tight_layout()
ComboFig.savefig('Combo_Teff_variance', dpi = 1200)
#%% Code to plot variance in results based on initial conditions
ControlFig = plt.figure()
Control = plt.subplot(111)
Control.plot(ts, CD4s_Control_Data_m_30pctEffector[1,:], 'g--', label = '30% Teff')
Control.plot(ts, CD4s_Control_Data_m_20pctEffector[1,:], 'g', label = '20% Teff')
Control.plot(ts, CD4s_Control_Data_m_10pctEffector[1,:], 'g:', label = '10% Teff')
Control.fill_between(ts, CD4s_Control_Data_m_30pctEffector[1,:],CD4s_Control_Data_m_10pctEffector[1,:],color = 'green', alpha = 0.2)
Control.plot(ts, CD4s_Control_Data_m_30pctEffector[2,:], 'r--', label = '30% Teff')
Control.plot(ts, CD4s_Control_Data_m_20pctEffector[2,:], 'r', label = '20% Teff')
Control.plot(ts, CD4s_Control_Data_m_10pctEffector[2,:], 'r:', label = '10% Teff')
Control.fill_between(ts, CD4s_Control_Data_m_30pctEffector[2,:],CD4s_Control_Data_m_10pctEffector[2,:],color = 'red', alpha = 0.2)
Control.set_xlabel('Time (Hours)', size = 22)
Control.set_ylabel('Fraction of CD4 T cells', size = 22)
Control.set_xticks(ticks = [0,50,100,150])
Control.set_xticklabels(labels = [0,50,100,150], fontsize = 18)
Control.set_yticks(ticks = [0,0.2,0.4,0.6,0.8,1.0])
Control.set_yticklabels(labels = [0,0.2,0.4,0.6,0.8,1.0], fontsize = 18)
Control.set_ylim(0,1.0)
Control.set_xlim(0,150)
ControlFig.tight_layout()
ControlFig.savefig('Control_Teff_variance', dpi = 1200)