#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:25:21 EDT 2019

@author: Bryan Medina
"""
from nwb_plots_functions import *
# READ ME ################################
'''
This file
- Gets the different values for t1, t2, ...,  t5, beta1, beta2, ..., beta5 for each trial
- Compares then all against each other.
'''
##########################################

## CHANGE ME #############################################################
# Data directory
DIRECTORY   = '/Users/bjm/Documents/CMU/Research/data/'
TRIAL_DATA  = '/Users/bjm/Documents/CMU/Research/data/trial_data/'
TRIAL_PLOTS = '/Users/bjm/Documents/CMU/Research/data/plots/trials/'
MOUSE_ID    = '424448'
##########################################################################


# Get file from directory
spikes_nwb_file = os.path.join(DIRECTORY, 'mouse' + MOUSE_ID + '.spikes.nwb')
nwb = h5.File(spikes_nwb_file, 'r')

probe_names = nwb['processing']

# I should probably make a trial data structure, then fill that data struct with with a trial number, the t set and the beta set, the configuration of that trial


# Changes depending on the trial.
start     = 0 #in second
end       = 2000 #in seconds

# time stamps ( this never changes )
# This is SPECIFICALLY for the 'drifting_gratings_2' stimulus
timestamps  = nwb['stimulus']['presentation']['drifting_gratings_2']['timestamps'].value
stim_orient = nwb['stimulus']['presentation']['drifting_gratings_2']['data'].value

PLOTTING = True
## For every region,
for probe_name in probe_names:

    print(probe_name)
    # File to get data from
    probe_filename = DIRECTORY + MOUSE_ID + "_" + probe_name

    try:
        with open(probe_filename, 'rb') as f:
            probe = pickle.load(f)

    except FileNotFoundError:

        for probe_name in probe_names:
            saveProbeData(MOUSE_ID, probe_name, nwb)

    ## For EVERY trial,
    for trial_number in range(len(timestamps)):
        print("Trial number %d" % trial_number)
        # Check if we have this file

        try: 
            trial_file = TRIAL_DATA + "/" +  MOUSE_ID + "/" + probe_name + "/tr_" + str(trial_number)
            with open(trial_file+")", 'rb') as t:
                tr = pickle.load(t)

        except FileNotFoundError:
            trial = timestamps[trial_number]
            freq  = stim_orient[trial_number][1]
            angle = stim_orient[trial_number][3]

            # Checking for 'nans'
            if not (str(freq) == "nan") or not (str(angle) == "nan"):
                freq  = int(freq)
                angle = int(angle)
                
                config = str(freq) + "_" + str(angle)
                
                ## go through every cell in that region,
                ## find out how that cell is behaving IN THE TRIAL'S TIME FRAME,
                ## and save that activity to a vector...
                
                ## do that for every trial... essentially make PSTHs for every trial...
                curr_trial = np.zeros((len(bins), 1))
                
                for cell in probe.getCellList():
                    spikes = nwb['processing'][probe_name]['UnitTimes'][str(cell)]['times'].value
                    stimulus_spikes = binarySearch(spikes, trial, 0, len(spikes)-1)
                    
                    if not (type(stimulus_spikes) == type(-1)):
                        stimulus_spikes = (stimulus_spikes - trial[0])
                        stimulus_spikes *= 1000
                        
                        for stim_spike in stimulus_spikes:
                            curr_trial[insertToBin(stim_spike, end)] += 1

                ########################
                tr        = Trial()
                tr.number = trial_number
                tr.config = config
                tr.spikes = curr_trial
                # tr.t
                # tr.beta

                z = fromFreqList(curr_trial)
                curr_trial,b,c = plt.hist(z, bins)
                plt.clf()

                curr_trial /= 0.001*len(probe.getCellList())

                tr.spikes = curr_trial
                tr.lsq    = LSQUnivariateSpline(bins[0:len(bins)-1], curr_trial, knots)
                #tr.lsq    = UnivariateSpline(bins[0:len(bins)-1], curr_trial)

                # Find the 
                #######################
                with open(trial_file, 'wb') as t:
                    pickle.dump(tr, t)

        if(PLOTTING):
            plt.xlim(-2, 500)
            plt.ylim(0, 50)
            plt.ylabel('Spikes/second')
            plt.xlabel('Bins')
            plt.title("Mouse: " + str(MOUSE_ID) + " | " + probe_name + " trial: " + str(tr.number) + " | " + tr.config)
            plt.bar(bins[0:len(bins)-1], tr.spikes, alpha=0.8, color='blue')
            plt.plot(xs, tr.lsq(xs), color='red', alpha=0.4)
            plt.show()
            #plt.savefig(TRIAL_PLOTS + MOUSE_ID + "/" + probe_name + "/" + "tr_"+str(trial_number))
            plt.clf()
            
