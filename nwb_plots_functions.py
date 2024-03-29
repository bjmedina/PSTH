#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:25:21 EDT 2019

@author: Bryan Medina 
"""
###### Imports ########
from random import choices
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.vectors import FloatVector
from scipy.interpolate import LSQUnivariateSpline
from scipy.interpolate import CubicSpline
from scipy.interpolate import UnivariateSpline

import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import sys
########################

# Setting up packages for rpy2 use
package_name = 'cpm'

if rpackages.isinstalled(package_name):
    have_package = True
    print("R package %s already installed" % package_name)
else:
    have_pakcage = False

if not have_package:
    utils = rpackages.importr('utils')
    utils.chooseCRANmirror(ind=1)

    utils.install_packages(package_name)
    print("installed R package: %s" % package_name)

cpm = rpackages.importr(package_name)
##################################

### This is for the 'drifting gratings' stimulus
# All possible temporal frequencies for the stimulus
temp_freqs   = [1, 2, 4, 8, 15]

# All possible orientations of stimulus (angles)
orientations = [i*45 for i in range(8)]

# Knots for spline (selected by eye)
knots    = [30, 50, 52, 55, 65, 70, 75, 80,  83, 100, 150, 200, 250, 300, 325, 375, 400]
tr_knots = [50, 110, 160, 200, 250, 300, 350, 400, 450]

# Number of times timulus is presented
num_trials = 600

# Conversion frmo ms to s
msToSec    = 1000 # 1000 ms in 1 sec

# For future plotting
xs = np.linspace(0,600,3000)

# Start and end of trials
start = 0
end   = 2000

# Bin width
width = 1

# Actual bins (for later use)
bins  = np.linspace(start, end, int( (end - start)/width + 1 ))

# Probe to region mapping
mapping = {'probeA': 'AM',
        'probeB': 'PM',
        'probeC': 'V1',
        'probeD': 'LM',
        'probeE': 'AL',
        'probeF': 'RL'}

colors = ['k', '#9400D3', 'b', 'g', '#FF7F00', 'r']

###

class Probe:

    # Max firing rate and the time it occurs
    max_frate  = 0
    max_ftime  = 0

    # Second highest firing rate
    max_frate2 = 0
    max_ftime2 = 0

    # Min firing rate and the time it occurs
    min_frate  = 0
    min_ftime  = 0

    # Second lowest 
    min_frate2 = 0
    min_ftime2 = 0

    # Average firing rate that is converged to as t -> 500 ms 
    converge   = 0

    # Change point (before the first peak)
    change_pt  = 0
    chg_time   = 0 

    # Average firing rate
    avg_frate  = 0

    # Standard deviation of the firing rates
    std        = 0

    # LSQUnivariate function 
    lsq        = " "
    
    def __init__(self, nwb, name):
        '''
        Description
        -----------
        Constructor
        
        Input(s)
        --------
        'nwb': h5py._hl.files.File. 'spikes.nwb' dataset. 
        'probe': string. name of probe.
        
        Output(s)
        ---------
        New 'Probe' object
        '''
        
        self.__cells = getProbeCells(nwb, name)
        self.name    = mapping[name]
        
    def getCell(self, cell_number):
        '''
        Description
        -----------
        Method returns dictionary of cell "at index 'cell_number'"
        
        Input(s)
        --------
        'cell_number': int. key of a corresponding cells
        
        Output(s)
        ---------
        Dictionary of cell 'cell_number'
        '''
        
        return self.__cells[cell_number]

    def getCellList(self):
        '''
        Description
        -----------
        Method returns dictionary of cells
        
        Output(s)
        ---------
        Dictionary of cell 'cell_number'
        '''

        return self.__cells.keys()

    def __repr__(self):
        '''
        Description
        -----------
        Method replaces default '__str__' with one that prints out average spiking rate, 2 maximum and 2 minimum firing rates, and the time in which they occur. 

        Output(s)
        ---------
        String to print.
        '''

        return "%s\t Avg: %3.2f Std: %3.2f | Max: %3.2f @ %d | Max2: %3.2f @ %d | Min: %3.2f @ %d | Min2: %3.2f @ %d | Converges to %3.2f | Change: %3.2f @ %d" % (self.name, self.avg_frate, self.std, self.max_frate, self.max_ftime, self.max_frate2, self.max_ftime2, self.min_frate, self.min_ftime, self.min_frate2, self.min_ftime2, self.converge, self.change_pt, self.chg_time)

    def __str__(self):
        '''
        Description
        -----------
        Method replaces default '__repr__' with one that's great for LaTeX-table making.

        Output(s)
        ---------
        String to print.
        '''

        return "%s & %3.2f & %3.2f & (%3.2f, %d) & (%3.2f, %d) & (%3.2f, %d) & (%3.2f, %d) & %3.2f & (%3.2f, %d)\\\\" % (self.name, self.avg_frate, self.std, self.max_frate, self.max_ftime, self.max_frate2, self.max_ftime2, self.min_frate, self.min_ftime, self.min_frate2, self.min_ftime2, self.converge, self.change_pt, self.chg_time)

def getProbeCells(nwb, probe):
    '''
    Description
    -----------
    'GetProbeCells' gets dataset and returns all cells, for a given probe, that are in the Visual Cortex.
    
    Input(s)
    --------
    'nwb': h5py._hl.files.File. 'spikes.nwb' dataset. 
    'probe': string. name of probe.
    
    Output(s)
    ---------
    'v_cells': dict. Dictionary that all cells that are in V.
    '''
    
    # Get all cells with activity in V
    cells   = nwb['processing'][probe]['unit_list'].value
    v_cells = {} 
    
    for cell in cells:
        region = nwb['processing'][probe]['UnitTimes'][str(cell)]['ccf_structure'].value.decode('utf-8')
        
        if region[0] == 'V' or region[0] == 'v':
            v_cells[cell] = Cell()
            
    return v_cells


class Cell:

    max_frate = 0
    max_ftime = 0
    avg_frate = 0
    std       = 0
    name      = " " 
    lsq       = " "
    # Change point (before the first peak)
    change_pt  = 0
    chg_time   = 0 
    
    def __init__(self):
        '''
        Description
        -----------
        Constructor
        
        Output(s)
        ---------
        New 'Cell' object
        '''
        self.__table    = makeTable()

    def getSpikes(self, config):
        '''
        Description
        -----------
        Method returns table for given cell
        
        Input(s)
        --------
        'config': string. key of dictionary.
        
        Output(s)
        ---------
        table at certain config
        '''
        return self.__table[config]

    def addSpike(self, config, spike, end):
        '''
        Description
        -----------
        Method adds 1 to spike counts
        
        Input(s)
        --------
        'config': string. key of dictionary.
        'spike' : time of spike in seconds
        'end'   : end of trial 
        
        Output(s)
        ---------
        table at certain config
        '''
        # Find out index spike needs to be in.
        bn = insertToBin(spike, end)
        
        # Add one to ongoing count.
        self.__table[config][bn] += 1

    def __str__(self):
        '''
        Description
        -----------
        Method replaces default '__str__' with one that prints out average spiking rate, 2 maximum and 2 minimum firing rates, and the time in which they occur. 

        Output(s)
        ---------
        String to print.
        '''
        return  "Max: %3.2f\t Avg: %3.2f\t Std: %3.2f" % (self.max_frate, self.avg_frate, self.std)

class Trial:
    # The trial number
    number = -1

    # The configuration
    config = ""

    # Should be five values for each of these
    t      = [None]*5
    beta   = [None]*5

    # Need to make this go from start to end ... This will hold the PSTH.
    spikes = np.zeros((len(bins), 1))

    lsq    = []

    def __add__(self, other_trial):
        '''
        Description
        -----------
        Method overrides '+' operator so that you can add two Trial objects

        Input(s)
        --------
        'other_trial': Trial. Another trial object
        
        Output(s)
        ---------
        sum of two trials (adds spiking histogram)
        '''
        pass
        

def makeTable():
    '''
    Description
    -----------
    'makeTable' creates a dictionary to keep track of time bins for each possible orientation of stimulus. One for each cell.
    
    Output(s)
    ---------
    'table': dict. Dictionary that contains orientation combination as key and all cells that are in V.
    '''

    bins  = np.linspace(start, end, int( (end - start)/width + 1 )) 
    
    # In this table, each key is a different configuration of the stimulus
    # and each row corresponds to spikes in a time bin.
    table = {}

    for freq in temp_freqs:
        
        for angle in orientations:

            config        = str(freq) + "_" + str(angle) 
            table[config] = np.zeros((len(bins), 1))

    
    return table

def binarySearch(spikes, interval, start, end):
    '''
    Description
    -----------
    'binarySearch' will find the index of a spike in a certain interval. Once it finds the index of a spike in the interval, it will try to find all the spikes that are in that interval. Essentially a modified take on the classic binary search algorithm.
    
    Input(s)
    --------
    'spikes'  : list. list of all spikes of a given neuron.
    'interval': list. current time interval of stimulus (usually about 2 seconds).
    'start'   : int. beginning
    'end'     : int. end


    Output(s)
    ---------
    list. Returns list of spikes in a given interval (first spike is found using binary search, the rest with the 'spikesInInterval' method.    
    '''
    
    if end >= 1:
    
        mid_point = midpoint(start, end)
        
        # If our spike is inside the interval, let's return the index
        if inside(spikes[mid_point], interval):
            return spikesInInterval(spikes, interval, mid_point)

        # If our spike is greater than (or less than) the interval, let's adjust checking bounds
        elif spikes[mid_point] > interval[1]:

            next_midpoint = midpoint(start, mid_point-1)

            # If this is true, then we're going to hit a recursion error...
            # We don't want that to happen.
            if mid_point == next_midpoint:
                return -1
            
            return binarySearch(spikes, interval, start, mid_point-1)

        elif spikes[mid_point] < interval[0]:
            
            next_midpoint = midpoint(mid_point+1, end)
            
            # If this is true, then we're going to hit a recursion error...
            # We don't want this.            
            if mid_point == next_midpoint:
                return -1
            
            return binarySearch(spikes, interval, mid_point+1, end)

    else:

        return -1


def spikesInInterval(spikes, interval, known):
    '''
    Description
    -----------
    'spikesInInterval' will find all spikes in a certain interval based on the index of one found in the interval.  
    
    Input(s)
    --------
    'spikes'  : list. list of all spikes of a given neuron.
    'interval': list. current time interval of stimulus (usually about 2 seconds).
    'known'   : int. Index in 'spikes' of a known spike in the interval.

    Output(s)
    ---------
    'spike_set': set. indices of all spikes in the interval. This is converted to a list when returned.
    '''

    # Index of known spike
    i         = known

    # Boolean variables we'll be using to determine if we're either checking 'above' or 'below' the known value.
    # 'DOWN' is true because we'll start by checking below the known spike
    DOWN      = True
    UP        = False

    # Set of spikes. We'll be using a set because 1) sets can't have duplicates and 2) checking for duplicates can be done in constant O(1) time. 
    spike_set = set()
  
    # We don't want to check out of bounds of the spikes list.
    while i > -1 and i < len(spikes):

        if inside(spikes[i], interval) and DOWN:
            spike_set.add(spikes[i])
            i = i - 1
            
        elif not inside(spikes[i], interval) and DOWN:
            i    = known + 1
            UP   = True
            DOWN = False
            
        elif inside(spikes[i], interval) and UP:
            spike_set.add(spikes[i])
            i = i + 1
            
        elif not inside(spikes[i], interval) and UP:
            break

    # Convert set to list, then return.
    return np.array(list(spike_set))
    

def inside(spike, interval):
    '''
    Description
    -----------
    'inside' will determine if a spike is in an interval.
 
    Input(s)
    --------
    'spikes'  : list. list of all spikes of a given neuron.
    'interval': list. current time interval of stimulus (usually about 2 seconds).

    Output(s)
    --------
    boolean. True if spike is in interval. False otherwise.
    '''
    
    return spike >= interval[0] and spike <= interval[1]


def midpoint(start_rate, end_rate):
    '''
    Description
    -----------
    'midpoint' will calculate midpoint between two points
 
    Input(s)
    --------
    'start_rate'  : int. beginning
    'end_rate'    : int. end

    Output(s)
    --------
    int. midpoint between 'start_rate' and 'end_rate'
    '''
    
    return int(start_rate + (end_rate - start_rate)/2)

                           
def insertToBin(spiketime, end):
    '''
    Description
    -----------
    'insertToBin' will bin that a spiketime belongs in
 
    Input(s)
    --------
    'spiketime'  : int. spike time in ms
    'end'    : int. end of trial

    Output(s)
    --------
    int. idx. Index that the spiketime belongs to
    '''    
    
    idx = int( (spiketime - (spiketime % width)) / width )
    
    if( idx > end ): 
        #print("spiketime " + str(spiketime) + "\tidx " + str(idx))
        idx = end

    return idx 



def saveProbeData(MOUSE_ID, probe_name, nwb):
    '''
    Description
    -----------
    'saveProbeData' save the data, using pandas, of a certain mouse given a certain probe.
 
    Input(s)
    --------
    'MOUSE_ID'  : int. ID of mouse we'll be looking at
    'probe_name': string. name of probe
    'nwb'       : h5py._hl.files.File. Dataset

    Output(s)
    --------
    None.
    '''
        
    # Changes depending on the trial.
    start     = 0 #in second
    end       = 2000 #in seconds
    
    # time stamps ( this never changes )
    # This is SPECIFICALLY for the 'drifting_gratings_2' stimulus
    timestamps  = nwb['stimulus']['presentation']['drifting_gratings_2']['timestamps'].value
    stim_orient = nwb['stimulus']['presentation']['drifting_gratings_2']['data'].value
    
    
    ## Adding spikes
    # Get all cells that are in V for every probe
    #print(probe_name)
    probe = Probe(nwb, probe_name)
    
    # Going to want to save this information later.
    filename = MOUSE_ID + "_" + probe_name
    
    # ...get every cell. Then...
    cells = probe.getCellList()
    
    # ... for every cell...
    for cell in cells:
                
        # (Getting current cell)
        curr_cell = probe.getCell(cell)
        
        # ...get the current cells spiking activity.
        spikes = nwb['processing'][probe_name]['UnitTimes'][str(cell)]['times'].value
        
        # For every occurrence of each kind of stimulus
        for i in range(len(timestamps)):
            
            # Extract interval of stimulus, temporal frequency of stimulus, and angle of stimulus.
            trial = timestamps[i]
            freq  = stim_orient[i][1]
            angle = stim_orient[i][3]
            
            # Checking for 'nans'
            if not (str(freq) == "nan") or not (str(angle) == "nan"):
                
                freq  = int(freq)
                angle = int(angle)
                
                # Convert freq and angle to something that can be used as an index. 
                config = str(freq) + "_" + str(angle) 
                
                # Search for all spikes that are in this time frame. 
                stimulus_spikes = binarySearch(spikes, trial, 0, len(spikes)-1)
                
                
                if not (type(stimulus_spikes) == type(-1)):
                    # questionable but should do the trick (to get everything between 0 and 2000 ms)
                    stimulus_spikes = (stimulus_spikes - trial[0])
                    
                    stimulus_spikes *= 1000
                    
                    # For all the spikes you just found, add them to the their respective bin.
                    for stim_spike in stimulus_spikes:
                        curr_cell.addSpike(config, stim_spike, end)
                        
    print("Saving to " + filename)
    
    with open(filename, 'wb') as f:
        pickle.dump(probe, f)
        

def fromFreqList(x):
    '''
    Description
    -----------
    'fromFreqList' converts frequency list to a list of repitions based on index. This is usefull for histograms.

    Example
    -------
    fromFreqList([2,1,4,2]) => [0,0,1,2,2,2,2,3,3]
 
    Input(s)
    --------
    'x': list of ints. 

    Output(s)
    --------
    'z': list of ints.
    '''    
    z = []
    for i in range(len(x)):
        y = [ i for ii in range(int(x[i])) ]
        for num in y:
            z.append(num)

    return z

def robj_to_dict(robj):
    '''
    Description
    -----------
    'robj_to_dict' converts an R object to a python dictionary
 
    Input(s)
    --------
    'robj': R object 

    Output(s)
    --------
    dictionary.

    Source
    ------
    https://medium.com/bigdatarepublic/contextual-changepoint-detection-with-python-and-r-using-rpy2-fa7d86259ba9
    
    '''    
    return dict(zip(robj.names, map(list, robj)))
