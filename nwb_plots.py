#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:25:21 EDT 2019

@author: Bryan Medina
"""
from nwb_plots_functions import *
# READ ME ################################
# This file plots
# - (1) PSTHs for every cell (averaged across all trials) as well as a smoothed curve
# - (2) PSTHs for every probe (averaged across all trials and all cells) as well as a smoothed curve
# - (3) Smoothed curve for every probe
##########################################

## CHANGE ME #############################################################
# Data directory
DIRECTORY = '/Users/bjm/Documents/CMU/Research/data'
SUMMARY_PLOTS_DIRECTORY = '/Users/bjm/Documents/CMU/Research/data/plots/'
VAR_DIREC = '/Users/bjm/Documents/CMU/Research/data/plots/variations/'
MOUSE_ID = '421338'
##########################################################################


# Get file from directory
spikes_nwb_file = os.path.join(DIRECTORY, 'mouse' + MOUSE_ID + '.spikes.nwb')
nwb = h5.File(spikes_nwb_file, 'r')

probe_names = nwb['processing']

# Allows plotting (takes more time)
PLOTTING = True

# Print Descriptions
DESCRIPTIONS = True

# Turn this on if it's your first time running this code.
ALL_PLOTS = True

if(ALL_PLOTS):
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
            
            for probe_name in probe_names:
                saveProbeData(MOUSE_ID, probe_name, nwb)
            
            print("Run again")
            sys.exit(1)
    
        # Summary of all activity across all cells in a probe.
        x = np.zeros((len(bins), 1))

        # Plotting (1) #####################
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
                
            # Convert cell spiking data to a format 'plt.hist' will like
            z = fromFreqList(curr_cell)
            curr_cell,b,c = plt.hist(z, bins)
            plt.clf()
            
            # Normalize
            curr_cell /= num_trials*0.001
            
            # Get some information on the cell such as max firing rate, avg, std, and name
            ################# Finding peaks and valleys #######################
            probe.getCell(cell).max_frate = max(curr_cell[0:500])
            probe.getCell(cell).max_ftime = np.where(curr_cell[0:500] == probe.getCell(cell).max_frate)[0][0]
            probe.getCell(cell).avg_frate = np.mean(curr_cell[0:500])
            probe.getCell(cell).std       = np.std(curr_cell[0:500])
            probe.getCell(cell).name      = cell
            
            # Also get the associated firing rate curve for the cell
            lsq = LSQUnivariateSpline(bins[0:len(bins)-1], curr_cell, knots)
            probe.getCell(cell).lsq = lsq

            cpm_result = cpm.detectChangePoint(FloatVector(lsq(curr_cell[0:probe.getCell(cell).max_ftime])), cpmType='Student', ARL0=1000)
            cpm_result = robj_to_dict(cpm_result)
            
            probe.getCell(cell).change_pt = lsq(cpm_result['changePoint'][0])
            probe.getCell(cell).chg_time  = cpm_result['changePoint'][0]
            ####################################################################
            
            if(DESCRIPTIONS):
                print("Cell " + str(cell) + " : " + str(probe.getCell(cell))) 
                
            # Plotting
            if(PLOTTING):
                # Plotting normalized cell activity
                cell_filename  = MOUSE_ID + "_cell" + str(cell)
                plt.axvline(x=probe.getCell(cell).chg_time, alpha=0.5, linestyle='--', color='magenta')
                plt.ylim(0, 75)
                plt.xlim(-20, 520)
                plt.ylabel('Spikes/second')
                plt.xlabel('Bins')
                plt.title("Mouse: " + str(MOUSE_ID) + " / " +  probe_name + " in "+ probe.name + ". Cell: " + str(cell))
                plt.plot(xs, lsq(xs), color = 'magenta', alpha=0.9) 
                plt.bar(b[0:len(b)-1], curr_cell)
                plt.savefig(CELL_PLOTS_DIRECTORY + cell_filename + ".png")
                plt.clf()    
            # End Plotting (1) ####################
        
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
        # First we find the first peak and the time it occurs at.
        probe.max_frate  = max(x[0:500]) 
        probe.max_ftime  = np.where(x[0:500] == probe.max_frate)[0][0]
        
        # Now first valley
        probe.min_frate  = min(x[0:probe.max_ftime]) 
        probe.min_ftime  = np.where(x[0:probe.max_ftime] == probe.min_frate)[0][0] 
        
        # Now second peak
        probe.max_frate2 = max(x[200:300]) 
        probe.max_ftime2 = np.where(x[200:300] == probe.max_frate2)[0][0] + 200
        
        # Last valley
        probe.min_frate2 = min(x[probe.max_ftime:probe.max_ftime2])
        probe.min_ftime2 = np.where(x[probe.max_ftime:probe.max_ftime2] == probe.min_frate2)[0][0] + probe.max_ftime
        
        # The value it converges towards the end.
        probe.converge   = min(x[probe.max_ftime2:500])
        
        # Average firing rate + standard deviation
        probe.avg_frate  = np.mean(x[0:500])
        probe.std        = np.std(x[0:500])

        # Smoothed Function
        lsq = LSQUnivariateSpline(bins[0:len(bins)-1], x, knots)
        probe.lsq = lsq
        
        # Get the change point here 
        cpm_result = cpm.detectChangePoint(FloatVector(lsq(xs[probe.min_ftime-5:probe.max_ftime+1])), cpmType='Student', ARL0=1000)
        cpm_result = robj_to_dict(cpm_result)
        
        # Set chnage point and change point time
        probe.change_pt = lsq(cpm_result['changePoint'][0]+probe.min_ftime-5)
        probe.chg_time  = cpm_result['changePoint'][0]+probe.min_ftime-5
        ###################################################################
        
        if(DESCRIPTIONS):
            print(str(probe))
            
        # Plotting (2) ###############################################
        if(PLOTTING):
            # Plotting
            plt.axvline(x=probe.chg_time, color='red', linestyle='--', alpha=0.7)
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
        # End Plotting (2) ###########################################

        

# Plotting (3) ###############################################
# Here, we'll plot all curves for every region for a given mouse.
probes = []

# First, lets order the probe in terms of the time in which the max firing rate occurs
for probe_name in probe_names:
    
    probe_filename = MOUSE_ID + "_" + probe_name

    with open(probe_filename, 'rb') as f:
        # Plotting all curves for every region for a given mouse.
        probe = pickle.load(f)

    probes.append(probe)

probes.sort(key=lambda x: x.max_ftime)

# Finally, we can plot
for i in range(0, len(probes)):

    probe = probes[i]
    plt.ylabel('Firing Rate (Spikes/second)')
    plt.xlabel('Bins (ms)')
    plt.ylim(0, 12)
    plt.xlim(-20, 500)
    plt.title("Mouse: " + str(MOUSE_ID) + " | Average Firing Rates")
    plt.plot(xs, probe.lsq(xs), label = probe.name, color=colors[i])

plt.legend()
plt.savefig(SUMMARY_PLOTS_DIRECTORY + str(MOUSE_ID) + ".png")
plt.clf()
# End Plotting (3) ###########################################

