#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 22:49:56 2019

@author: L. Hauser (HAUZERO)
"""

### LOAD LIBRARIES ###
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.signal import savgol_filter
from sewar.full_ref import sam as sam

### FUNCTION TO FIND THE BEST FIT BETWEEN FIELD SPECTRA AND PROSAIL SOIL SPECTRA BASED ON PSOIL*RSOIL SIMULATIONS ####
def spectrasolver(x,soil_spectrum1, soil_spectrum2, rho_soil):
    psoil, rsoil    = x 
    sim_soil        = rsoil*(psoil*soil_spectrum1+(1-psoil)*soil_spectrum2)
    wavel           = range(400,2501,1)
    rho_soil1 = np.transpose(np.array([wavel, rho_soil]))
    sim_soil1 = np.transpose(np.array([wavel, sim_soil]))
    samerror        = sam(rho_soil1, sim_soil1)
    
    return samerror

#### LOCATIONS ANALYZED ###
#[5, 14, 23, 57, 61, 65, 66, 68, 70]
locations = [5, 14, 23, 61, 65, 66]
saver = np.zeros((len(locations)+1,2101))

### STARTING OPTMIZATION WITH THESE PSOIL AND RSOIL SETTINGS AND RANGES ###
#PSOIL (soil moisture factor: 0 wet, 1 dry)
psoil = 0.7
rsoil = 1
bdns    = ((0.00,1.00), (0,1))
    
##### 2 EXTREME SOIL SPECTRA FROM PROSAIL WET VERSUS DRY ##########
soilspectra      = pd.read_csv('/media/leon/FREECOM HDD/Data/Portugal/soilspectraps.txt', delimiter = ' ', header=None)
soil_spectrum1   = soilspectra[0]
soil_spectrum2   = soilspectra[1]
#rho_soil = rsoil*(psoil*soil_spectrum1+(1-psoil)*soil_spectrum2

######## IN LOOP GOING OVER ALL LOCATIONS AND MEASUREMENTS TO SELECT THE ADEQUATE MEASUREMENTS ####
for ilc, loca in enumerate(locations):
    print 'location: ' + str(locations[ilc])
    path = '/media/leon/FREECOM HDD/Data/Portugal/Portugal2019_Fielddata/Location ' + str(loca) + '/Soil Spectra/'
    files = glob.glob(path + '*.sed')
    
    for ix, filex in enumerate(files):
        spec = pd.read_csv(filex, skiprows=26, sep='\t')
        if ix == 0:
            savedspec   = np.array(spec['Wvl'])
        plt.plot(spec['Wvl'], spec['Reflect. %'])
        plt.show()
        bla = raw_input("Save? y/n")
        if bla == 'y':
            savedspec = np.vstack((savedspec, spec['Reflect. %']))
            print 'saved!'
    saver[0]        = savedspec[0,50:2151]
    saver[ilc+1]    = np.mean(savedspec[1:5,50:2151], axis=0)

### PLOTTING AVERAGE SPECTRA PER LOCATION ######
for ig in range(len(saver)):
    fig = plt.figure()
    plt.title('location: ' + str(locations[ig-1]))
    plt.plot(saver[0], saver[ig])
    plt.show

#### AVERAGE OVER ALL LOCATIONS + SMOOTHING #########
rho_soil =  (np.mean(saver[1:,:], axis=0)/100.)
rho_soil =  savgol_filter(rho_soil,19,1)
plt.plot(saver[0], rho_soil)

#### CALLING THE OPTIMIZATION FUNCTION #######
x = psoil, rsoil
res = minimize(spectrasolver,x,args=(soil_spectrum1, soil_spectrum2, rho_soil), method='nelder-mead',bounds=bdns)
psoil,rsoil = res.x
sim_soil        = rsoil*(psoil*soil_spectrum1+(1-psoil)*soil_spectrum2)
plt.plot(saver[0], sim_soil, label='Optimized Sim')
plt.plot(saver[0], rho_soil, label='Actual Spectrum')
plt.plot(saver[0], soilspectra[1], label='Wet Sim')
plt.plot(saver[0], soilspectra[0], label='Dry Sim')
plt.legend()