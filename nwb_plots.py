#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:25:21 EDT 2019

@author: Bryan Medina
"""
###### Imports ########
from nwb_plots_functions import *
from scipy.interpolate import LSQUnivariateSpline

import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import sys
########################


## CHANGE ME ####################################################################################
# Data directory
DIRECTORY = '/Users/bjm/Documents/CMU/Research/data'

VAR_DIREC = '/Users/bjm/Documents/CMU/Research/data/plots/variations/'

MOUSE_ID = '405751'
#################################################################################################

# Get file from directory
spikes_nwb_file = os.path.join(DIRECTORY, 'mouse' + MOUSE_ID + '.spikes.nwb')
nwb = h5.File(spikes_nwb_file, 'r')

probe_names = nwb['processing']

# Length of stimulus and desired bin width
start      = 0 # in ms
end        = 2000 #in ms
width      = 1 # ms
msToSec    = 1000 # 1000 ms in 1 sec

# Actual bins (for later use)
bins       = np.linspace(start, end, int( (end - start)/width + 1 ))

# Number of times timulus is presented
num_trials = 630 
    
# All possible orientations of stimulus (angles and temporal frequencies) [drifting gratings]
orientations = [i*45 for i in range(8)]
temp_freqs   = [1, 2, 4, 8, 15]

# For future plotting
xs = np.linspace(0,2000,10000)

# Knots for spline (selected by eye)
knots = [50, 70, 100, 150, 200, 250, 300, 325, 375, 400]

# Allows plotting (takes more time)
PLOTTING = True

# Print Descriptions
DESCRIPTIONS = True

# The median, 90th percentile neuron, and 10th percentile neuron.
median_neuron = []
neuron_90th   = []
neuron_10th   = []

for probe_name in probe_names:
    # File to get data from.
    probe_filename = MOUSE_ID + "_" + probe_name
    print(probe_filename)
    
    # plot directories
    
    ## CHANGE ME ####################################################################################
    PROBE_PLOTS_DIRECTORY = '/Users/bjm/Documents/CMU/Research/data/plots/probes/'
    CELL_PLOTS_DIRECTORY  = '/Users/bjm/Documents/CMU/Research/data/plots/cells/' + probe_name + '/'
    #################################################################################################

    ## Find probe to override
    try:
        with open(probe_filename, 'rb') as f:
            probe = pickle.load(f)
    ## If probe file doesn't exist, then we'll have to make that file from scratch        
    except FileNotFoundError:
        saveProbeData(MOUSE_ID, probe_name, nwb)
        print("Run again")
        sys.exit(1)
        

    # Summary of all activity across all cells in a probe.
    x = np.zeros((len(bins), 1))
    
    # Getting all data for a given cell
    for cell in probe.getCellList():
        # current cell spiking data
        curr_cell = np.zeros((len(bins), 1))
        for freq in temp_freqs:
            for angle in orientations:
                config = str(freq) + "_" + str(angle)
                curr_cell += probe.getCell(cell).getSpikes(config)
                # Plot curr cell
                x += probe.getCell(cell).getSpikes(config)

        
        z = fromFreqList(curr_cell)
        curr_cell,b,c = plt.hist(z, bins)
        plt.clf()
        curr_cell /= num_trials*0.001
        probe.getCell(cell).max_frate = max(curr_cell[0:500])
        probe.getCell(cell).avg_frate = np.mean(curr_cell[0:500])
        probe.getCell(cell).std       = np.std(curr_cell[0:500])
        probe.getCell(cell).name      = cell


        if(DESCRIPTIONS):
            print(probe.getCell(cell)) 
        
        lsq = LSQUnivariateSpline(bins[0:len(bins)-1], curr_cell, knots[1:-1])
        probe.getCell(cell).lsq = lsq
        
        # Plotting
        if(PLOTTING):
            # Plotting normalized cell activity
            cell_filename  = MOUSE_ID + "_cell" + str(cell)
            print(cell_filename)
            plt.ylim(0, 75)
            plt.xlim(-20, 520)
            plt.ylabel('Spikes/second')
            plt.xlabel('Bins')
            plt.title("Mouse: " + str(MOUSE_ID) + " / " +  probe_name + " in "+ probe.name + ". Cell: " + str(cell))
            plt.plot(xs, lsq(xs), color = 'magenta', alpha=0.76) 
            plt.bar(b[0:len(b)-1], curr_cell)
            plt.savefig(CELL_PLOTS_DIRECTORY + cell_filename + ".png")
            plt.clf()
    

    # Plotting normalized probe activity
    z = fromFreqList(x)
    
    x,b,c = plt.hist(z, bins)
    plt.clf()
    ###
    
    ### Normalization
    # also divide by number of neurons in that particular region
    x /= num_trials*(0.001)*len(probe.getCellList())

    # Need to find the two maxes and two mins

    ################# Finding peaks and valleys #######################
    # First we find the first peak
    probe.max_frate  = max(x[0:500])
    probe.max_ftime  = np.where(x[0:500] == probe.max_frate)[0][0]

    # Now first valley
    probe.min_frate  = min(x[0:probe.max_ftime])
    probe.min_ftime  = np.where(x[0:probe.max_ftime] == probe.min_frate)[0][0]

    # Now second peak
    probe.max_frate2 = max(x[200:300])
    probe.max_ftime2 = np.where(x[200:300] == probe.max_frate2)[0][0]

    # Last valley
    probe.min_frate2 = min(x[probe.max_ftime:probe.max_ftime2])
    print("x[probe.max_ftime:probe.max_ftime2]")
    probe.min_ftime2 = np.where(x[probe.max_ftime:probe.max_ftime2] == probe.min_frate2)[0][0]
    
    probe.avg_frate  = np.mean(x[0:500])
    probe.std        = np.std(x[0:500])
    ###################################################################

    if(DESCRIPTIONS):
        print(probe)

    
    # Smoothing 
    lsq = LSQUnivariateSpline(bins[0:len(bins)-1], x, knots[1:-1])
    probe.lsq = lsq

    if(PLOTTING):
        # Plotting
        plt.ylim(0, 12)
        plt.xlim(-20, 500)
        plt.ylabel('Spikes/second')
        plt.xlabel('Bins')
        plt.title("Mouse: " + str(MOUSE_ID) + " / " +  probe_name + " in "+ probe.name)
        plt.plot(xs, lsq(xs), color = 'red') 
        plt.bar(b[0:len(b)-1], x, alpha=0.8)
        plt.savefig(PROBE_PLOTS_DIRECTORY + probe_filename + ".png")
        
        plt.clf()

    with open(probe_filename, 'wb') as f:
        pickle.dump(probe, f)
'''
TODO: 
MORT1dSMOOTH
15.2.3 (leading upto as well)
BARS (poisson version)
419
'''
