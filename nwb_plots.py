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
BIN_WIDTH = 100 # in s

# open NWB files using h5py library

###### UPDATE PATH #################################
DIRECTORY = '/Users/bjm/Documents/CMU/Research/data'

# UPDATE MOUSE ID #
MOUSE_ID = '424448'

start     = 0 #in second
end       = 3*60*60 #in seconds

# Get file from directory
spikes_nwb_file = os.path.join(DIRECTORY, 'mouse' + MOUSE_ID + '.spikes.nwb')
nwb = h5.File(spikes_nwb_file, 'r')

# names of probes
names = nwb['processing'].keys()

# dictionary of probes
probes = {}

for probe_name in names:
    # Get all cells that are in V for every probe
    probes[probe_name] = Probe(nwb, probe_name, BIN_WIDTH, start, end)


# time stamps ( this never changes )
# This is SPECIFICALLY for the 'drifting_gratings_2' stimulus
timestamps  = nwb['stimulus']['presentation']['drifting_gratings_2']['timestamps'].value
stim_orient = nwb['stimulus']['presentation']['drifting_gratings_2']['data'].value

###

## Adding spikes to time bin
for probe_name in names:

    print(probe_name)
    
    cells = probes[probe_name].getCellList()

    for cell in cells:
        print("\t" + str(cell))
        curr_cell = probes[probe_name].getCell(cell)

        spikes = nwb['processing'][probe_name]['UnitTimes'][str(cell)]['times'].value

        for i in range(len(timestamps)):
            
            interval = timestamps[i]
            freq  = stim_orient[i][1]
            angle = stim_orient[i][3]
            
            if not (str(freq) == "nan") or not (str(angle) == "nan"):
                freq  = int(freq)
                angle = int(angle)
                
            orient = str(freq) + "_" + str(angle)

            stimulus_spikes = binarySearch(spikes, interval, 0, len(spikes)-1)

            if not (stimulus_spikes == -1):

                for stim_spike in stimulus_spikes:
                    curr_cell.addSpike(orient, insertToBin(stim_spike, BIN_WIDTH))


        
            
'''
TODO: (plotting)
'''
