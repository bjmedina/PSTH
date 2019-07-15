#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 11:50:02 EDT 2019

@author: Bryan Medina
"""
###### Imports ########
from nwb_plots_functions import *
########################

###### UPDATE PATH #################################
DIRECTORY = '/Users/bjm/Documents/CMU/Research/data'
VAR_DIREC = '/Users/bjm/Documents/CMU/Research/data/plots/variations/'
MOUSE_ID  = '424448'
MICE_ID = ['424448', '421338', '405751']
####################################################

# Get file from directory
spikes_nwb_file = os.path.join(DIRECTORY, 'mouse' + MOUSE_ID + '.spikes.nwb')
nwb = h5.File(spikes_nwb_file, 'r')

probe_names = nwb['processing']

# keeps track of max firing rate for each cell in 
probe_fr = {}

alphas = {'424448': 1,
          '421338': 0.6,
          '405751': 0.4}

# firing rate filename
filename = MOUSE_ID + '_probes_fr'
PLOT_ALL = False

# We want to check to see if we have this data
try:
    with open(filename, 'rb') as f:
        probe_fr = pickle.load(f)
        
except:
    # only keep track of maximal firing rates...
    probe_fr = {}
    
    for probe_name in probe_names:
        # Getting all data for a given cell
        # File to get data from.
        probe_filename = MOUSE_ID + "_" + probe_name
        print(probe_filename)
        
        try:
            with open(probe_filename, 'rb') as f:
                # Plotting all curves for every region for a given mouse.
                probe = pickle.load(f)
                
        except FileNotFoundError:
            saveProbeData(MOUSE_ID, probe_name, nwb)
            print("Run again nwb_plots with plotting off")
            sys.exit(1)
            
        probe_fr[probe_name] = []
        
        for cell in probe.getCellList():
            # Get max, add it here...
            probe_fr[probe_name].append(probe.getCell(cell).max_frate)


# Plot everything
for probe_name in probe_names:
    # Plot variability of every region
    if(PLOT_ALL):
        # Plotting how variable neuron can be
        for probe_name in probe_names:

            plt.title("Mouse: " + str(MOUSE_ID) + " / " + probe_name + " Variation")
            plt.ylim(0, 14)
            plt.xlabel("Maximal Firing Rate (Spikes/Sec)")
            plt.ylabel("Number of Neurons")
            plt.hist(probe_fr[probe_name], bins = 100)
            plt.savefig(VAR_DIREC + MOUSE_ID + probe_name +  "_variations.png")
            plt.clf()

for MOUSE in MICE_ID:
    
    variability = []

    filename = MOUSE + '_probes_fr'

    with open(filename, 'rb') as f:
        probe_fr = pickle.load(f)

    for probe_name in probe_names:
        variability.extend(probe_fr[probe_name])
                
    plt.hist(variability, bins = 100, label="Mouse: " + MOUSE, alpha=alphas[MOUSE], edgecolor='black')

plt.legend()
plt.title("Variability of Maximal Firing Rates")
plt.xlabel("Maximal Firing Rate (Spikes/Sec)")
plt.ylabel("Number of Neurons")
plt.savefig(VAR_DIREC + "variations.png")

# Save the probe_fr file.
with open(filename, 'wb') as f:
    pickle.dump(probe_fr, f)


