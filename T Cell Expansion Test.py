# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:14:27 2020

@author: davem
"""

from scipy import integrate
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
#These import modules that I use a lot. You will pretty much ALWAYS use numpy and matploylib.pyplot
#The other ones are because I'm going to be doing a lot of differential equation solving, so I
#imported integrate. I probably don't need interpolate.

def ThreePopModel(ys, t0, params):
    #define the function that is our model to produce derivative for integration
    [p_R, p_E, T_R, T_E, k_RE] = params
#    import the paramteres
    N0 = ys[0]
    R0 = ys[1]
    E0 = ys[2]
#    import the starting values and assign them to variable that make more sense
    
    dydt = np.zeros(3)
#    create an array to hold our derivative values which will later be the output of the function
    
    dNdt = 0.0 - N0/p_R - N0/p_E
    dRdt = N0/p_R + R0/T_R
    dEdt = N0/T_E + (E0 - k_RE*R0)/T_E
#    solve for the actual values of our derivatives
    
    dydt[0] = dNdt
    dydt[1] = dRdt
    dydt[2] = dEdt
#    store those values in the array that we created so that we can easily export them from this function
    return dydt
    

def ExpansionSim(f_R0, f_E0, runTime = 6.0*24, cells_tot_0 = 100000.0, TGF = False, Rapa = False, species = 'mouse',stepsize = 1):
#     initialization of expansion simulator. cells_tot_0 is the total number of cells initially
#    which will usually be  100,000. f_r0 is the initial
#    fraction of cells that are regulatory T cells, f_e0 is the initial fraction
#     of cells that are effector T cells. The rest of the cells will be naive T cells
#   including TGF as true or rapa as true will shift certain constants and which differential
#   equations are used. Setting the species to 'human' will eventually do the same. 
    N = cells_tot_0*(1.0-f_R0-f_E0)
    R = cells_tot_0*f_R0
    E = cells_tot_0*f_E0
#    N is the number of naive T cells, R is the number of regulatory T cells, E is the number of 
#    effector T cells
    if species == 'mouse':
#        if loop to specify what species we are dealing with
        if TGF == False:
#            if loop to determine if there is TGF-beta present. 
            p_R = float('inf')
            p_E = 1.5*24 #days*hours to get half life in hours
#            P_R and P_E are parameters that determine the differentiation rate of naive t cells into 
#            the different phenotypes. values should be ~1/2 lives in hours
        elif TGF == True:
#           if there is TGF 
            p_R = 1.5*24#days*hours to get half life in hours
            p_E = 1.5*24 #days*hours to get half life in hours
#            P_R and P_E are parameters that determine the differentiation rate of naive t cells into 
#            the different phenotypes. values should be ~1/2 lives in hours
        else: #here to check error if you forgot to specify TGF presence
            print 'No TGF parameter specified.'
        if Rapa == False:
            T_R = 2*24.0 #days*hours to get dbling time in hours
            T_E = 1.25*24.0 #days*hours to get dbling time in hours
#            T_R and T_E are parameters that determine the doubling time of the respective
#            phenotypes
        elif Rapa == True:
            T_R = 2*24.0
            T_E = 3*24.0
        else:
            print 'No Rapa parameter specified.'
        k_RE = 1.0
#        k_RE is a parameter that determines how many effector T cells a regulatory
#        T cell suppresses on average. Starting with 1, will dig into literature
    elif species == 'human':
#        if loop to change parameters if the species is human
        if TGF == False:
#            if loop to determine if there is TGF-beta present. 
            p_R = 0
            p_E = 1.5*24 #days*hours to get half life in hours
#            P_R and P_E are parameters that determine the differentiation rate of naive t cells into 
#            the different phenotypes. values should be ~1/2 lives in hours
        elif TGF == True:
#           if there is TGF 
            p_R = 1.5*24#days*hours to get half life in hours
            p_E = 1.5*24 #days*hours to get half life in hours
#            P_R and P_E are parameters that determine the differentiation rate of naive t cells into 
#            the different phenotypes. values should be ~1/2 lives in hours
        else: #here to check error if you forgot to specify TGF presence
            print 'No TGF parameter specified.'
        if Rapa == False:
            T_R = 3*24.0 #days*hours to get dbling time in hours
            T_E = 2*24.0 #days*hours to get dbling time in hours
#            T_R and T_E are parameters that determine the doubling time of the respective
#            phenotypes
        elif Rapa == True:
            T_R = 3*24.0
            T_E = 4*24.0
        else:
            print 'No Rapa parameter specified.'
        k_RE = 1.0
#        k_RE is a parameter that determines how many effector T cells a regulatory
#        T cell suppresses on average. Starting with 1, will dig into literature
    else:
        print 'Species not supported'
    params = [p_R, p_E, T_R, T_E, k_RE]
    #package parameters to import into the diff EQ model. I might make the parameter generator a
    #separate piece later on. 
    
    y0s = [N, R, E]
    ts = np.arange(0, runTime + 1, stepsize)  
#    package our starting y values (y0s) and our t values (ts) to put into ODE solver   
    
    yOut = integrate.odeint(ThreePopModel, y0s, ts, args=(params,))
    
    Ns = yOut[:,0]
    Rs = yOut[:,1]
    Es = yOut[:,2]
    
    return [Ns, Rs, Es, ts, params]
    
    
def FractionCalcs (Ns, Rs, Es)