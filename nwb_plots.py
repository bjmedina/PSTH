#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:25:21 EDT 2019

@author: Bryan Medina
"""
###### Imports ########
from nwb_plots_functions import *
from scipy.interpolate import CubicSpline
from scipy.ndimage.filters import gaussian_filter1d

import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
########################

MOUSE_ID   = '424448'
probe_name = 'probeA'
start      = 0 # in ms
end        = 2000 #in ms
width      = 1 # ms
bins       = np.linspace(start, end, int( (end - start)/width + 1 ))
num_trials = 630 # You should just get this from the data...
    
# All possible orientations of stimulus (angles and temporal frequencies)
orientations = [i*45 for i in range(8)]
temp_freqs   = [1, 2, 4, 8, 15]

# File to get data from.
filename = MOUSE_ID + "_" + probe_name

# plot directory
PLOTS = '/Users/bjm/Documents/CMU/Research/data/plots/'


try:
    with open(filename, 'rb') as f:
        probe = pickle.load(f)
        
except FileNotFoundError:
    saveData(MOUSE_ID)
    print("Run again")

# Summary of all activity across all cells in a probe.
x = np.zeros((len(bins), 1))

# Getting all data
for cell in probe.getCellList():
    for freq in temp_freqs:
        for angle in orientations:
            config = str(freq) + "_" + str(angle)
            x = x + probe.getCell(cell).getSpikes(config)

### Improve this
z = []
for i in range(len(x)):
    y = [ i for ii in range(int(x[i])) ]
    for num in y:
        z.append(num)

x,b,c = plt.hist(z, bins)
plt.clf()
###

### Normalization
# also divide by number of neurons in that particular region
x /= num_trials*(0.0001)*len(probe.getCellList()) #Should I divide by number of cells?

# Smoothing (with cubic spline)
cs = CubicSpline(bins[0:len(bins)-1], x)    
plt.plot(bins, cs(bins)) 


# Plotting
plt.title(filename + ": " + probe.name)
plt.bar(b[0:len(b)-1], x)
plt.savefig(PLOTS + filename + ".png")


'''smoothing
-----------
'''
