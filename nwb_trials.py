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
TMAX_DIREC  = '/Users/bjm/Documents/CMU/Research/data/tmax/'
MOUSE_ID    = '421338'
##########################################################################


# Get file from directory
spikes_nwb_file = os.path.join(DIRECTORY, 'mouse' + MOUSE_ID + '.spikes.nwb')
nwb = h5.File(spikes_nwb_file, 'r')

probe_names = nwb['processing']

# Whether or not we want to calculate the confidence intervals
CONF_INTERVAL = True
BOOTSTRAPS = 500

trials = []
t_max1 = []
t_max2 = []

# Changes depending on the trial.
start     = 0 #in second
end       = 2000 #in seconds

# time stamps ( this never changes )
# This is SPECIFICALLY for the 'drifting_gratings_2' stimulus
timestamps  = nwb['stimulus']['presentation']['drifting_gratings_2']['timestamps'].value
stim_orient = nwb['stimulus']['presentation']['drifting_gratings_2']['data'].value

PLOTTING = False
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
            with open(trial_file, 'rb') as t:
                tr = pickle.load(t)
                trials.append(tr)

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

                trials.append(tr)
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
            


    if(CONF_INTERVAL):
        # Calculating the confidence intervals
        fname = TMAX_DIREC + MOUSE_ID + "/" + probe_name + "/" + MOUSE_ID + "_tmax_"

        try:
            with open(fname + "1", 'rb') as f:
                t_max1 = pickle.load(f)

            with open(fname + "2", 'rb') as f:
                t_max2 = pickle.load(f)

        except FileNotFoundError:
            # We're doing 500 bootstraps
            for i in range(0, BOOTSTRAPS):
                print("BOOTSTRAP %d" % i)
                # g is going to be our random sample, size 600, of the 600 trials
                g = choices(trials, k = len(trials))
                
                sample_spikes = np.zeros((len(g[0].spikes),))
                lsq = np.zeros((len(g[0].lsq(xs)), 1))
                
                # Now we need to construct our curves based on these 600 samples
                for sample in g:
                    # Need to add all spikes together
                    ## To do this, we have to *essentially* do an element wise addition
                    for j in range(0, len(sample.spikes)):
                        sample_spikes[j] += sample.spikes[j]

                # Recompute tmax_1 and tmax_2
                ## We have to normalize sample_spikes by number of trials
                sample_spikes /= len(g)
                peak  = max(sample_spikes[0:500])
                tmax_1 = np.where(sample_spikes[0:500] == peak)[0][0]
                
                peak2 = max(sample_spikes[200:300])
                tmax_2 = np.where(sample_spikes[200:300] == peak2)[0][0] + 200
                
                if(PLOTTING):
                    print("Peak 1: %d @ %d" % (peak, tmax_1))
                    print("Peak 2: %d @ %d" % (peak, tmax_2))
                    plt.ylim(0, 10)
                    plt.xlim(-2, 500)
                    plt.bar(bins[:-1], sample_spikes, alpha=0.8, color='blue')
                    plt.axvline(x=tmax_1,color='red', linestyle='--')
                    plt.axvline(x=tmax_2,color='red', linestyle='--')
                    
                    plt.show()
                    plt.clf()
            
                # Save those two into two separate vectors
                t_max1.append(tmax_1)
                t_max2.append(tmax_2)
            

    
    # clear the slate for the next probe         
    trials = []
    
    with open(fname + "1", 'wb') as f:
        pickle.dump(t_max1, f)

        
    with open(fname + "2", 'wb') as f:
        pickle.dump(t_max2, f)

    t_max1 = []
    t_max2 = []

    
fname = TMAX_DIREC + MOUSE_ID + "/" + probe_name + "/" + MOUSE_ID + "_tmax_"
