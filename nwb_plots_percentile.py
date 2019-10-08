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
PERC_PLOTS_DIRECTORY = '/Users/bjm/Documents/CMU/Research/data/plots/percentile/'
MOUSE_ID = '424448'
####################################################

# Get file from directory
spikes_nwb_file = os.path.join(DIRECTORY, 'mouse' + MOUSE_ID + '.spikes.nwb')
nwb = h5.File(spikes_nwb_file, 'r')

probes = nwb['processing']
probe_names = [name for name in probes.keys()]

# save all curves for all regions
mid = {}
top = {}
bot = {}

# Used for plotting
rows = 3
cols = 2

for probe_name in probe_names:
    # Calculate median neuron, and also 90th and 10th percentile neuron
    median_n = []
    top_ten  = []
    bot_ten  = []
    
    probe_filename = MOUSE_ID + "_" + probe_name

    with open(probe_filename, 'rb') as f:
        probe = pickle.load(f)
        
    for xval in xs:
        
        rates = []
        
        for cell in probe.getCellList():
            rates.append(probe.getCell(cell).lsq(xval))

        # Sort this list...
        rates.sort()
        
        median_n.append(np.median(rates))
        top_ten.append(np.percentile(rates, 75))
        bot_ten.append(np.percentile(rates, 25))
     
    # save the curves
    mid[probe_name] = LSQUnivariateSpline(xs, median_n, knots[1:-1])
    top[probe_name] = LSQUnivariateSpline(xs, top_ten, knots[1:-1])
    bot[probe_name] = LSQUnivariateSpline(xs, bot_ten, knots)


# Plotting median, 75th percentile, and 25th percentile neuron

# Do multiple plots on one figure
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(10, 10))
fig.tight_layout(pad=0.1, w_pad=0.1, h_pad=0.1)
fig.suptitle("Mouse %s Neural Activity" % (MOUSE_ID))
fig.text(0.5, 0.04, 'Bins (ms)', ha='center')
fig.text(0.04, 0.5, 'Firing Rate (Spike/sec)', va='center', rotation='vertical')
i = 0

for row in range(0, rows):
    for col in range(0, cols):

        probe_name = probe_names[i]
        probe_filename = MOUSE_ID + "_" + probe_name
        
        with open(probe_filename, 'rb') as f:
            probe = pickle.load(f)

        box = axes[row,col].get_position()
        move = 0.08
        move2 = 0.033
        move3 = 0.053

        if(row == 0):
            if(col == 0):
                axes[row,col].set_position([move+box.x0+box.x0/5, box.y0, box.width * 0.8 , box.height * 0.8])
            else:
                axes[row,col].set_position([move+box.x0-box.x0/7, box.y0, box.width * 0.8 , box.height * 0.8])
        elif(row == 1):
            if(col == 0):
                axes[row,col].set_position([move+box.x0+box.x0/5, box.y0+move2, box.width * 0.8 , box.height * 0.8])
            else:
                axes[row,col].set_position([move+box.x0-box.x0/7, box.y0+move2, box.width * 0.8 , box.height * 0.8])
        elif(row == 2):
            if(col == 0):
                axes[row,col].set_position([move+box.x0+box.x0/5, box.y0+move3, box.width * 0.8 , box.height * 0.8])
            else:
                axes[row,col].set_position([move+box.x0-box.x0/7, box.y0+move3, box.width * 0.8 , box.height * 0.8])

                
        axes[row, col].set_ylim([0, 13])
        axes[row, col].set_xlim([-20, 500])
            
        axes[row, col].set_title(probe.name)
        
        axes[row, col].plot(xs, top[probe_name](xs), label = "75th Percentile")
        
        axes[row, col].plot(xs, mid[probe_name](xs), label = "Median Neuron")
    
        axes[row, col].plot(xs, bot[probe_name](xs), label = "25th Percentile")

        if(row == 0 and col == cols - 1):
            axes[row, col].legend()
        
        # Next probe
        i = i+1

plt.savefig(PERC_PLOTS_DIRECTORY + str(MOUSE_ID) + "_percentile.png")
plt.clf()

