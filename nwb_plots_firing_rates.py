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
MICE_ID = ['424448', '421338', '405751']
MOUSE_ID  = '421338'

####################################################

# Get file from directory
spikes_nwb_file = os.path.join(DIRECTORY, 'mouse' + MOUSE_ID + '.spikes.nwb')
nwb = h5.File(spikes_nwb_file, 'r')

probe_names = nwb['processing']

# keeps track of max firing rate for each cell in 
probe_fr = {}

colors = {'424448':'red',
          '421338':'green',
          '405751':'blue'}

# firing rate filename
filename = MOUSE_ID + '_probes_fr'
PLOT_ALL = True

rows = 2
cols = 2

# Ideally, you should do this for every mouse.

# We want to check to see if we have this data
try:
    with open(filename+"_", 'rb') as f:
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
            plt.hist(probe_fr[probe_name], bins = 100, edgecolor='black')
            plt.savefig(VAR_DIREC + MOUSE_ID + probe_name +  "_variations.png")
            plt.clf()


# Plotting multiple summary plots in one plot.
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8,8))
fig.suptitle("Variation in Maximal Firing Rates")
fig.text(0.5, 0.04, 'Maximal Firing Rate (Spikes/sec)', ha='center')
fig.text(0.04, 0.5, 'Number of Neurons', va='center', rotation='vertical')

variability = []
curves = {}
i = 0

# Plotting 4 plots in one figure.
for row in range(0, rows):
    for col in range(0, cols):

        if( not (row + 1 == rows and col + 1 == cols) ):
            MOUSE = MICE_ID[i]
            filename = MOUSE + '_probes_fr'
            
            with open(filename, 'rb') as f:
                probe_fr = pickle.load(f)
                
            for probe_name in probe_names:
                variability.extend(probe_fr[probe_name])
                
            axes[row, col].set_ylim([0, 90])
            axes[row, col].set_xlim([0, 100])
            axes[row, col].set_title("Mouse %s" % (MOUSE))
            ys, bins, c = axes[row, col].hist(variability, bins = 100,color=colors[MOUSE], edgecolor='black', alpha=0.7)    
            curves[MOUSE] = [LSQUnivariateSpline(bins[0:len(bins)-1], ys, [10, 30, 55, 70, 100]), bins[0:len(bins)-1]]
            i = i+1
            variability = []

        else:
            axes[row, col].set_ylim([0, 90])
            axes[row, col].set_xlim([0, 100])
            axes[row, col].set_title("All Variations")

            for ID in MICE_ID:
                axes[row, col].plot(curves[ID][1], curves[ID][0](curves[ID][1]), label=ID, color=colors[ID], alpha=0.7)
            axes[row, col].legend()

plt.savefig(VAR_DIREC + "firing_rate_variations.png")

# Save the probe_fr file.
with open(filename, 'wb') as f:
    pickle.dump(probe_fr, f)


