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
DIRECTORY = '/Users/bjm/Documents/CMU/Research/data'
MOUSE_ID = '424448'
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

# For every region, 
# For EVERY trial,
# go through every cell in that region, find out how that cell is behaving IN THE TRIAL'S TIME FRAME, and save that activity to a vector... do that for every trial... essentially make PSTHs for every trial... should only be 600...


