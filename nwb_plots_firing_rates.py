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

MOUSE_ID = '421338'

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

        ## Scrap this idea... only keep track of maximal firing rates...
        probe_fr[probe_name] = {}
        probe_fr[probe_name]["high"]   = 0
        probe_fr[probe_name]["medium"] = 0
        probe_fr[probe_name]["low"]    = 0
        
        try:
            with open(probe_filename, 'rb') as f:
                probe = pickle.load(f)
                
        except FileNotFoundError:
            saveProbeData(MOUSE_ID, probe_name, nwb)
            print("Run again")
            sys.exit(1)

        # Scrap this
        high = probe.avg_frate + 2*probe.std
        low  = probe.avg_frate - 2*probe.std
            
        # Getting all data for a given cell
        for cell in probe.getCellList():
            # current cell spiking data
            curr_cell = np.zeros((len(bins), 1))
            for freq in temp_freqs:
                for angle in orientations:
                    config = str(freq) + "_" + str(angle)
                    curr_cell += probe.getCell(cell).getSpikes(config)
                    
            z = fromFreqList(curr_cell)
            curr_cell,b,c = plt.hist(z, bins)
            plt.clf()
            curr_cell /= num_trials*0.001

            #print("HI: %f \t LO: %f" % (high, low))
            for spike_rate in curr_cell:
                if spike_rate >= high:
                    probe_fr[probe_name]["high"] += 1
                    #print("HI\t\t%f\t\tprobe_fr[probe_name][high]:%d" % (spike_rate, probe_fr[probe_name]["high"]))
                elif spike_rate <= low:
                    probe_fr[probe_name]["low"] += 1
                    #print("LOW\t\t%f\t\tprobe_fr[probe_name][low]:%d" % (spike_rate, probe_fr[probe_name]["low"]))
                else:
                    probe_fr[probe_name]["medium"]+=1
                    #print("MID\t\t%f\t\tprobe_fr[probe_name][mid]:%d" % (spike_rate, probe_fr[probe_name]["medium"]))

    with open(filename, 'wb') as f:
        pickle.dump(probe_fr, f)


for probe_name in probe_names:

    plt.title("Mouse: " + str(MOUSE_ID) + " / " + probe_name + " Variation")

    plt.bar(probe_fr[probe_name].keys(), probe_fr[probe_name].values())
    plt.savefig(VAR_DIREC + MOUSE_ID + probe_name +  "_variations.png")
    plt.clf()
