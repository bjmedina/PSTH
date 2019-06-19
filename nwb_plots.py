#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:25:21 EDT 2019

@author: Bryan Medina
"""
###### Imports ########
from nwb_plots_functions import *
from scipy.ndimage.filters import gaussian_filter1d

import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
import os
########################


#### UPDATE SIZE ####
BIN_WIDTH = 1000 # in s

# open NWB files using h5py library

###### UPDATE PATH #################################
DIRECTORY = '/Users/bjm/Documents/CMU/Research/data'

# UPDATE MOUSE ID #
MOUSE_ID = '424448'

start     = 0 #in second
end       = 3*60*60 #in seconds

bins      = [(start+i*BIN_WIDTH) for i in range(end-start+1)]
num_bins  = len(bins)

# Get file from directory
spikes_nwb_file = os.path.join(DIRECTORY, 'mouse' + MOUSE_ID + '.spikes.nwb')
nwb = h5.File(spikes_nwb_file, 'r')

# names of probes
names = nwb['processing'].keys()

# dictionary of probes
probes = {}

for probe_name in names:
    # Get all cells that are in V for every probe
    probes[probe_name] = Probe(nwb, probe_name)


# time stamps ( this never changes )
# This is SPECIFICALLY for the 'drifting_gratings_2' stimulus
timestamps  = nwb['stimulus']['presentation']['drifting_gratings_2']['timestamps'].value
stim_orient = nwb['stimulus']['presentation']['drifting_gratings_2']['data'].value

###

## Adding spikes to time bin

# For every probe...
for probe_name in names:

    print(probe_name)

    # ...get every cell. Then...
    cells = probes[probe_name].getCellList()

    # ... for every cell...
    for cell in cells:
        
        print("\t" + str(cell))
        
        # (Getting current cell)
        curr_cell = probes[probe_name].getCell(cell)

        # ...get the current cells spiking activity.
        spikes = nwb['processing'][probe_name]['UnitTimes'][str(cell)]['times'].value

        # For every occurrence of each kind of stimulus
        for i in range(len(timestamps)):

            # Extract interval of stimulus, temporal frequency of stimulus, and angle of stimulus.
            interval = timestamps[i]
            freq  = stim_orient[i][1]
            angle = stim_orient[i][3]

            # Checking for 'nans'
            if not (str(freq) == "nan") or not (str(angle) == "nan"):
                freq  = int(freq)
                angle = int(angle)

            # Convert freq and angle to something that can be used as an index. 
            orient = str(freq) + "_" + str(angle)

            # Search for all spikes that are in this time frame. 
            stimulus_spikes = binarySearch(spikes, interval, 0, len(spikes)-1)

            if not (stimulus_spikes == -1):

                # For all the spikes you just found, add them to the their respective bin.
                for stim_spike in stimulus_spikes:
                    curr_cell.addSpike(orient, stim_spike)

freq = probes['probeA'].getCell(43).getSpikes("1_0")

'''
TODO: probe summer + cell summary (do cell summary first)
'''


'''
TODO: (plotting)
'''
temp_freqs   = [1, 2, 4, 8, 15]
orientations = [i*45 for i in range(8)]

for freq in temp_freqs:
    for angle in orientations:
        orient = str(freq) + "_" + str(angle)
        data = probes['probeA'].getCell(43).getSpikes(orient)
        plt.hist(data, bins)
