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

# Number of different stimulus combinations possible
# {1, 2, 4, 8, 15} X {0, 45, 90, 135, 180, 225, 270, 315} X 15 + 1 kind of grey image (30 of them)
DIFF_COMBS = 41

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


### things you can use

# gets spiking times for cell 43 in probeA ( loop throught probes and cell numbers )
spikes = nwb['processing']['probeA']['UnitTimes']['43']['times'].value

# time stamps ( this never changes )
timestamps = nwb['stimulus']['presentation']['drifting_gratings_2']['timestamps'].value


