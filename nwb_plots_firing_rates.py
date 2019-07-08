#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 11:50:02 EDT 2019

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

###### UPDATE PATH #################################
DIRECTORY = '/Users/bjm/Documents/CMU/Research/data'
VAR_DIREC = '/Users/bjm/Documents/CMU/Research/data/plots/variations/'

MOUSE_ID = '424448'

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
    
# All possible orientations of stimulus (angles and temporal frequencies)
orientations = [i*45 for i in range(8)]
temp_freqs   = [1, 2, 4, 8, 15]


probe_fr = {}

filename = MOUSE_ID + '_probes_fr'

try:
    with open(filename+"_", 'rb') as f:
        probe_fr = pickle.load(f)
except:

    for probe_name in probe_names:
        # File to get data from.
        probe_filename = MOUSE_ID + "_" + probe_name
        print(probe_filename)

        ## f only keep track of maximal firing rates...
        probe_fr[probe_name] = []
        
        try:
            with open(probe_filename, 'rb') as f:
                probe = pickle.load(f)
                
        except FileNotFoundError:
            saveProbeData(MOUSE_ID, probe_name, nwb)
            print("Run again")
            sys.exit(1)
        
        # Getting all data for a given cell
        for cell in probe.getCellList():
            # Get max, add it here...
            probe_fr[probe_name].append(probe.getCell(cell).max_frate)

    with open(filename, 'wb') as f:
        pickle.dump(probe_fr, f)


for probe_name in probe_names:

    plt.title("Mouse: " + str(MOUSE_ID) + " / " + probe_name + " Variation")

    plt.hist(probe_fr[probe_name])
    plt.savefig(VAR_DIREC + MOUSE_ID + probe_name +  "_variations.png")
    plt.clf()
