# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:14:27 2020

@author: davem
"""
#### V1.1 Update: I added in differences between CD4 and CD8 parameters for the mouse model. Still need to add the differences into the human model
#### V1.2 Update: I added in the parameters for mouse T cells according to the data that I took from flow. Calculations for parameters were done by fitting the
    #analytical solutions to various data points in Excel. See


from scipy import integrate
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
#These import modules that I use a lot. You will pretty much ALWAYS use numpy and matploylib.pyplot
#The other ones are because I'm going to be doing a lot of differential equation solving, so I
#imported integrate. I probably don't need interpolate.

#%% This is the function that defines the mathematical model we are using to model population dynamics
def ThreePopModelMouse(ys, t0, params):
    #define the function that is our model to produce derivative for integration
    [p_R, p_E, T_R, T_E8, T_E4, k_RE] = params
#    import the paramteres
    N40 = ys[0]
    R0 = ys[1]
    E80 = ys[2]
    E40 = ys[3]
#    import the starting values and assign them to variable that make more sense
    
    dydt = np.zeros(4)
#    create an array to hold our derivative values which will later be the output of the function

#    solve for the actual values of our derivatives
    dNdt = 0.0 - N40/p_R - N40/p_E
    dRdt = N40/p_R + R0/T_R
    if E40 - k_RE*R0 >= 0:
        dE4dt = N40/p_E + (E40 - k_RE*R0)/T_E4
    else:
        dE4dt = N40/p_E
#        I threw in ^this^ if loop to account for the fact that Treg suppression shouldn't cause
#        Teffectors to decrease, but could potentially put it at a standstill. Maybe having the
#        retarding term be a divisor rather than subtractor makes more sense
    if E80 - k_RE*R0 >= 0:
        dE8dt = (E80 - k_RE*R0)/T_E8
    else:
        dE8dt = 0 

    
    dydt[0] = dNdt
    dydt[1] = dRdt
    dydt[2] = dE8dt
    dydt[3] = dE4dt
#    store those values in the array that we created so that we can easily export them from this function
    return dydt
   
#%% This is the function that defines the mathematical model we are using to model population dynamics.
#   The difference between the human and mouse model is that the human model has a delay built in for 
#   the expansion parameters (no expansion before 48 hours) but still has differentiation
def ThreePopModelHuman(ys, t0, params):
    #define the function that is our model to produce derivative for integration
    [p_R, p_E, T_R, T_E8, T_E4, k_RE] = params
#    import the paramteres
    N40 = ys[0]
    R0 = ys[1]
    E80 = ys[2]
    E40 = ys[3]
#    import the starting values and assign them to variable that make more sense
    
    dydt = np.zeros(4)
#    create an array to hold our derivative values which will later be the output of the function
    
#    solve for the actual values of our derivatives. In the human model I am just shifting the start of expansion to after 2 days.
    if t0 > 48:
        dNdt = 0.0 - N40/p_R - N40/p_E
        dRdt = N40/p_R + R0/T_R
        if E40 - k_RE*R0 >= 0:
            dE4dt = N40/p_E + (E40 - k_RE*R0)/T_E4
        else:
            dE4dt = N40/p_E
    #        I threw in ^this^ if loop to account for the fact that Treg suppression shouldn't cause
    #        Teffectors to decrease, but could potentially put it at a standstill. Maybe having the
    #        retarding term be a divisor rather than subtractor makes more sense
        if E80 - k_RE*R0 >= 0:
            dE8dt = (E80 - k_RE*R0)/T_E8
        else:
            dE8dt = 0 
#        I threw in ^this^ if loop to account for the fact that Treg suppression shouldn't cause
#        Teffectors to decrease, but could potentially put it at a standstill. Maybe having the
#        retarding term be a divisor rather than subtractor makes more sense
    else:
        dNdt = 0.0 - N40/p_R - N40/p_E
        dRdt = N40/p_R
        dE4dt = N40/p_E
        dE8dt = 0
        
#    dNdt = 0.0 - N40/p_R - N40/p_E
#    dRdt = N40/p_R + R0/T_R
#    if E40 - k_RE*R0 >= 0:
#        dE4dt = N40/p_E + (E40 - k_RE*R0)/T_E4
#    else:
#        dE4dt = N40/p_E
##        I threw in ^this^ if loop to account for the fact that Treg suppression shouldn't cause
##        Teffectors to decrease, but could potentially put it at a standstill. Maybe having the
##        retarding term be a divisor rather than subtractor makes more sense
#    if E80 - k_RE*R0 >= 0:
#        dE8dt = (E80 - k_RE*R0)/T_E8
#    else:
#        dE8dt = 0

    
    dydt[0] = dNdt
    dydt[1] = dRdt
    dydt[2] = dE8dt
    dydt[3] = dE4dt
#    store those values in the array that we created so that we can easily export them from this function
    return dydt   

#%% This is the main function I will be using that compiles and runs the simulation
def ExpansionSim(f_R0, f_E80, f_E40, runTime = 6.0*24, cells_tot_0 = 100000.0, TGF = False, Rapa = False, species = 'mouse',stepsize = 1):
#     initialization of expansion simulator. cells_tot_0 is the total number of cells initially
#    which will usually be  100,000. f_r0 is the initial
#    fraction of cells that are regulatory T cells, f_e40 is the initial fraction
#     of cells that are effector CD4 T cells. f_E80 is the intial fraction that are CD8 T cells. 
#     The rest of the cells will be naive T cells
#   including TGF as true or rapa as true will shift certain constants and which differential
#   equations are used. Setting the species to 'human' will eventually do the same. 
    N = cells_tot_0*(1.0-f_R0-f_E80-f_E40)
    R = cells_tot_0*f_R0
    E8 = cells_tot_0*f_E80
    E4 = cells_tot_0*f_E40
#    N is the number of naive T cells, R is the number of regulatory T cells, E is the number of 
#    effector T cells
#%% This section is for MOUSE parameters
    if species == 'mouse':
#        if loop to specify what species we are dealing with
        if TGF == False:
#            if loop to determine if there is TGF-beta present. 
            p_R = float('inf') #Naive T cell "half-life" to become regulatory CD4 in hours. Infinite in the absence of TGF-beta
            p_E4 = 15.63460135 #Naive T cell "half-life" to become effector CD4 in hours
#            P_R and P_E4 are parameters that determine the differentiation rate of naive t cells into 
#            the different phenotypes. values should be ~1/2 lives in hours
        elif TGF == True:
#           if there is TGF 
            p_R = 17.73696497 #Naive T cell "half-life" to become effector CD4 in hours
            p_E4 = 15.63460135 #days*hours to get half life in hours
#            P_R and P_E are parameters that determine the differentiation rate of naive t cells into 
#            the different phenotypes. values should be ~1/2 lives in hours
        else: #here to check error if you forgot to specify TGF presence
            print ('No TGF parameter specified.')
        if Rapa == False:
            T_R = 44.75272703 #dbling time in hours
            T_E4 = 49.59029057 #dbling time in hours
            T_E8 = 34.53550219# dbling time in hours
#            T_R and T_E are parameters that determine the doubling time of the respective
#            phenotypes
        elif Rapa == True:
            T_R = 44.75272703
            T_E4 = 61.23003942
            T_E8 = 51.54913153
        else:
            print ('No Rapa parameter specified.')
        k_RE = 0.00
#        k_RE is a parameter that determines how many effector T cells a regulatory
#        T cell suppresses on average. Starting with 1, will dig into literature
#%% This section is for HUMAN parameters
    elif species == 'human':
#        if loop to change parameters if the species is human
        if TGF == False:
#            if loop to determine if there is TGF-beta present. 
            p_R = float('inf')
            p_E4 = 15.63460135 #days*hours to get half life in hours
#            P_R and P_E are parameters that determine the differentiation rate of naive t cells into 
#            the different phenotypes. values should be ~1/2 lives in hours
        elif TGF == True:
#           if there is TGF 
            p_R = 17.73696497 #days*hours to get half life in hours
            p_E4 = 15.63460135 #days*hours to get half life in hours
#            P_R and P_E are parameters that determine the differentiation rate of naive t cells into 
#            the different phenotypes. values should be ~1/2 lives in hours
        else: #here to check error if you forgot to specify TGF presence
            print ('No TGF parameter specified.')
        if Rapa == False:
            # T_R = 31.75272703 #dbling time in hours
            # T_E4 = 55.59029057 #dbling time in hours
            # T_E8 = 60.53550219# dbling time in hours #days*hours to get dbling time in hours
            T_R = 44.75272703 #dbling time in hours
            T_E4 = 49.59029057 #dbling time in hours
            T_E8 = 34.53550219# dbling time in hours

#            T_R and T_E are parameters that determine the doubling time of the respective
#            phenotypes
        elif Rapa == True:
            # T_R = 50.75272703
            # T_E4 = 186.23003942
            # T_E8 = 197.
            T_R = 44.75272703
            T_E4 = 61.23003942
            T_E8 = 51.54913153
        else:
            print ('No Rapa parameter specified.')
        k_RE = 0.00
#        k_RE is a parameter that determines how many effector T cells a regulatory
#        T cell suppresses on average. Starting with 1, will dig into literature
    else:
        print ('Species not supported')
#%% This section packages the parameters and runs the appropriate model depending on whether or not it is human or mouse.
    
    params = [p_R, p_E4, T_R, T_E8, T_E4, k_RE]
    #package parameters to import into the diff EQ model. I might make the parameter generator a
    #separate piece later on. 
    
    y0s = [N, R, E8, E4]
    ts = np.arange(0, runTime + 1, stepsize)  
#    package our starting y values (y0s) and our t values (ts) to put into ODE solver   
    if species == 'mouse':
        yOut = integrate.odeint(ThreePopModelMouse, y0s, ts, args=(params,))
#    this outputs a matrix of size (ts,ys), where yOut[:,0] will give the value for this first
#    variable (in this case N) for all timesteps. The extra comma in args=(params,) converts it
#    into a tuple which is necessary to feed it to args
    elif species == 'human':
        yOut = integrate.odeint(ThreePopModelHuman, y0s, ts, args=(params,))
    else:
        print ('Species not supported')
    
    Ns = yOut[:,0]
    Rs = yOut[:,1]
    E8s = yOut[:,2]
    E4s = yOut[:,3]
#    reassign outputs to variables that make sense
#    Have the function return what we want.
    return [Ns, Rs, E8s, E4s, ts, params]
    
#%% This is a simple function to help with some of the conversion of data from the previous section into fractional values for graphing
def FractionCalcs(Ns, Rs, E4s, E8s):
#    function for calculating the fractional portion of each population at all time steps
    cell_totals = Ns + Rs + E4s + E8s
    f_E8 = E8s/cell_totals
    f_R = Rs/cell_totals
    f_N = Ns/cell_totals
    f_E4 = E4s/cell_totals
    
#    calculations to find total cell count and fractional populations
    
    return [cell_totals, f_N, f_E8, f_E4, f_R]
    
    
#%% To see how these work together uncomment (highlight the relevant text and hit Ctrl + 1) and run the following:

## This first set of following code will simulate the expansion of mouse T cells
## over the course of 7 days from a starting population of 100,000 cells, 1% of 
## which are Tregs, 16% are effector CD4, and 45.2% are CD8s.
#    
#[Ns, Rs, E8s, E4s, ts, params] = ExpansionSim(0.011, 0.452, 0.163)
#
#[cell_totals, f_N, f_E8, f_E4, f_R] = FractionCalcs(Ns, Rs, E4s, E8s)
#
#countsFig = plt.figure()
#counts = plt.subplot(111)
#counts.plot(ts, E8s, 'k', label = 'CD8 T Cells')
#counts.plot(ts, Ns,'b', label = 'Naive CD4 T Cells')
#counts.plot(ts, E4s, 'r', label = 'Effector CD4 T Cells')
#counts.plot(ts, Rs, 'g', label = 'Regulatory T Cells')
#counts.legend(loc = 'best')
#
#freqsFig = plt.figure()
#freqs = plt.subplot(111)
#freqs.plot(ts, f_E8, 'k', label = 'CD8 T Cells')
#freqs.plot(ts, f_N,'b', label = 'Naive CD4 T Cells')
#freqs.plot(ts, f_E4, 'r', label = 'Effector CD4 T Cells')
#freqs.plot(ts, f_R, 'g', label = 'Regulatory T Cells')
#freqs.legend(loc = 'best')
#
##I'm going to save data generated from this as an np.array which just makes it easier to save
## and re-import later.
#Control_Data = np.array([Ns, Rs, Es, f_N, f_R, f_E, ts])

#f_tot_CD4s = f_N + f_E4 + f_R
#f_N_CD4s = f_N/f_tot_CD4s
#f_R_CD4s = f_R/f_tot_CD4s
#f_E4_CD4s = f_E4/f_tot_CD4s
#
#CD4s_Control_Data = np.array([f_N_CD4s, f_R_CD4s, f_E4_CD4s, ts])

## This second set of code will simulate the the same expansion as before, but now there
## is TGF beta in the media, which I have accounted for only in its ability to enhance
## naive T cell conversion, but not in suppression of effector T cells. I also have a variable
## that accounts for Treg suppression of effector T cell expansion, but right now I have it 
## set to zero, and I might just leave it out for now.


#[Ns, Rs, E8s, E4s, ts, params] = ExpansionSim(0.011, 0.452, 0.163, TGF = True)
#
#[cell_totals, f_N, f_E8, f_E4, f_R] = FractionCalcs(Ns, Rs, E4s, E8s)
#
#countsFig = plt.figure()
#counts = plt.subplot(111)
#counts.plot(ts, E8s, 'k', label = 'CD8 T Cells')
#counts.plot(ts, Ns,'b', label = 'Naive CD4 T Cells')
#counts.plot(ts, E4s, 'r', label = 'Effector CD4 T Cells')
#counts.plot(ts, Rs, 'g', label = 'Regulatory T Cells')
#counts.legend(loc = 'best')
#
#freqsFig = plt.figure()
#freqs = plt.subplot(111)
#freqs.plot(ts, f_E8, 'k', label = 'CD8 T Cells')
#freqs.plot(ts, f_N,'b', label = 'Naive CD4 T Cells')
#freqs.plot(ts, f_E4, 'r', label = 'Effector CD4 T Cells')
#freqs.plot(ts, f_R, 'g', label = 'Regulatory T Cells')
#freqs.legend(loc = 'best')
#
#TGF_Data= np.array([Ns, Rs, Es, f_N, f_R, f_E, ts])

#f_tot_CD4s = f_N + f_E4 + f_R
#f_N_CD4s = f_N/f_tot_CD4s
#f_R_CD4s = f_R/f_tot_CD4s
#f_E4_CD4s = f_E4/f_tot_CD4s
#
#CD4s_TGF_Data = np.array([f_N_CD4s, f_R_CD4s, f_E4_CD4s, ts])

## This third set of code will simulate the the same expansion as before, but now there
## is only rapamycin, which effectively slows down the expansion of effector T cells, but not
## the regulatory T cells. However, there is nothing inducing regulatory T cells in this one,
## so the population dynamics of the Tregs should be the same as the first simulation.

#[Ns, Rs, E8s, E4s, ts, params] = ExpansionSim(0.011, 0.452, 0.163, Rapa = True)
#
#[cell_totals, f_N, f_E8, f_E4, f_R] = FractionCalcs(Ns, Rs, E4s, E8s)
#
#countsFig = plt.figure()
#counts = plt.subplot(111)
#counts.plot(ts, E8s, 'k', label = 'CD8 T Cells')
#counts.plot(ts, Ns,'b', label = 'Naive CD4 T Cells')
#counts.plot(ts, E4s, 'r', label = 'Effector CD4 T Cells')
#counts.plot(ts, Rs, 'g', label = 'Regulatory T Cells')
#counts.legend(loc = 'best')
#
#freqsFig = plt.figure()
#freqs = plt.subplot(111)
#freqs.plot(ts, f_E8, 'k', label = 'CD8 T Cells')
#freqs.plot(ts, f_N,'b', label = 'Naive CD4 T Cells')
#freqs.plot(ts, f_E4, 'r', label = 'Effector CD4 T Cells')
#freqs.plot(ts, f_R, 'g', label = 'Regulatory T Cells')
#freqs.legend(loc = 'best')
#
#Rapa_Data = np.array([Ns, Rs, Es, f_N, f_R, f_E, ts])

#f_tot_CD4s = f_N + f_E4 + f_R
#f_N_CD4s = f_N/f_tot_CD4s
#f_R_CD4s = f_R/f_tot_CD4s
#f_E4_CD4s = f_E4/f_tot_CD4s
#
#CD4s_Rapa_Data = np.array([f_N_CD4s, f_R_CD4s, f_E4_CD4s, ts])

## This Last set of code will simulate the expansion with both Rapa and TGF beta present
## and should take into account both suppression of effector T cell expansion due to rapa
## as well as phenotypic conversion to Tregs due to TGF.

#[Ns, Rs, E8s, E4s, ts, params] = ExpansionSim(0.011, 0.452, 0.163,TGF = True, Rapa = True)
#
#[cell_totals, f_N, f_E8, f_E4, f_R] = FractionCalcs(Ns, Rs, E4s, E8s)
#
#countsFig = plt.figure()
#counts = plt.subplot(111)
#counts.plot(ts, E8s, 'k', label = 'CD8 T Cells')
#counts.plot(ts, Ns,'b', label = 'Naive CD4 T Cells')
#counts.plot(ts, E4s, 'r', label = 'Effector CD4 T Cells')
#counts.plot(ts, Rs, 'g', label = 'Regulatory T Cells')
#counts.legend(loc = 'best')
#
#freqsFig = plt.figure()
#freqs = plt.subplot(111)
#freqs.plot(ts, f_E8, 'k', label = 'CD8 T Cells')
#freqs.plot(ts, f_N,'b', label = 'Naive CD4 T Cells')
#freqs.plot(ts, f_E4, 'r', label = 'Effector CD4 T Cells')
#freqs.plot(ts, f_R, 'g', label = 'Regulatory T Cells')
#freqs.legend(loc = 'best')
#
#Combo_Data = np.array([Ns, Rs, Es, f_N, f_R, f_E, ts])

#f_tot_CD4s = f_N + f_E4 + f_R
#f_N_CD4s = f_N/f_tot_CD4s
#f_R_CD4s = f_R/f_tot_CD4s
#f_E4_CD4s = f_E4/f_tot_CD4s
#
#CD4s_Combo_Data = np.array([f_N_CD4s, f_R_CD4s, f_E4_CD4s, ts])


#%% Feel free to mess around with the input values for the ExpansionSim.
#Break it in different ways and see what it does. It isn't very robust, and I didn't include
#a lot of error messages. 