#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:25:21 EDT 2019

@author: Bryan Medina 

"""
###### Imports ########
from scipy.ndimage.filters import gaussian_filter1d

import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
import os
########################

class Probe:

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
        Method returns dictionary of cell "at index 'cell_number'"
        
        Input(s)
        --------
        'cell_number': int. key of a corresponding cells
        
        Output(s)
        ---------
        Dictionary of cell 'cell_number'
        '''

        return self.__cells.keys()


class Cell:
    
    def __init__(self, start, end, bin_width, nwb, probe):
        '''
        Description
        -----------
        Constructor
        
        Input(s)
        --------
        'start'    : integer. start time of experiment in seconds
        'end'      : integer. end time of experiment in seconds
        'bin_width': integer. width of each time bin in seconds (should this all be in seconds?)
        'num_combs': integer. number of different combinations the stimuli can take on.
    
        
        Output(s)
        ---------
        New 'Cell' object
        '''
        self.__table = makeTable(start, end, bin_width, nwb, probe)
        
    
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

    ## Calculating number of bins needed based on the predefined bin width
    ########
    ##### Is this variable??????
    #######
    start     = 0 #in seconds
    end       = 3*60*60 #in seconds
    bin_width = 100
    
    for cell in cells:
        region = nwb['processing'][probe]['UnitTimes'][str(cell)]['ccf_structure'].value.decode('utf-8')
        
        if region[0] == 'V' or region[0] == 'v':
            v_cells[cell] = Cell(start, end, bin_width, nwb, probe)
            
    return v_cells


def makeTable(start, end, bin_width, nwb, probe):
    '''
    Description
    -----------
    'makeTable' creates a dictionary to keep track of time bins for each possible orientation of stimulus.

    Input(s)
    --------
    'start'    : integer. start time of experiment in seconds
    'end'      : integer. end time of experiment in seconds
    'bin_width': integer. width of each time bin in seconds (should this all be in seconds?)
    'nwb': h5py._hl.files.File. 'spikes.nwb' dataset. 
    'probe': string. name of probe.
    
    Output(s)
    ---------
    'table': dict. Dictionary that contains orientation combination as key and all cells that are in V.
    '''
    num_combs = 41

    # Get number of bins based on bin width, and start/end time of experiment
    num_bins = int((end-start)/bin_width+1)

    # Corresponds to times of experiment
    bins     = np.linspace(start, end, num_bins, dtype=int)
  
    # All possible temporal frequencies for the stimulus
    temp_freqs   = [1, 2, 4, 8, 15]
    
    # All possible orientations of stimulus (angles)
    orientations = [i*45 for i in range(8)]
    
    # In this table, each key is a different configuration of the stimulus
    # and each row corresponds to spikes in a time bin.
    table        = {}

    for freq in temp_freqs:
        
        for angle in orientations:

            config        = str(freq) + "_" + str(angle)
            table[config] = np.zeros((num_bins, 2))
            
            table[config][:, 0] = bins

    # Empty images
    table['NaN_NaN'] = np.zeros((num_bins, 2))
    table['NaN_NaN'][:, 0] = bins
     
    return table


def makeProbeTable(probes, start, end, bin_width, num_combs):
    '''
    Description
    -----------
    'makeProbeTable' creates a dictionary for each probe to keep track of time bins for each possible orientation of stimulus.

    Input(s)
    --------
    'Probe'    : list. list of names of all probes
    'start'    : integer. start time of experiment in seconds
    'end'      : integer. end time of experiment in seconds
    'bin_width': integer. width of each time bin in seconds (should this all be in seconds?)
    'num_combs': integer. number of different combinations the stimuli can take on.
    
    
    Output(s)
    ---------
    'tables': dict. Dictionary that contains probe cell name as key.
    '''
    
    tables = {}
    
    for probe in probes:
        
        tables[probe] = makeTable(start, end, bin_width, num_combs)

    return tables


def binarySearch(spikes, interval, start, end):
    '''
    Description
    -----------
    'binarySearch' will find the index of a spike in a certain interval. Once it finds the index of a spike in the interval, it will try to find all the spikes that are in that interval. Essentially a modified take on the classic binary search algorithm.
    
    Input(s)
    --------
    'spikes'  : list. list of all spikes of a given neuron.
    'interval': list. current time interval of stimulus (usually about 2 seconds).

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

            # If this is true, then we're going to hit a recursion error....
            # We don't want this.
            if mid_point == next_midpoint:
                return -1
            
            return binarySearch(spikes, interval, start, mid_point - 1)

        elif spikes[mid_point] < interval[0]:
            
            next_midpoint = midpoint(mid_point+1, end)
            
            # If this is true, then we're going to hit a recursion error....
            # We don't want this.            
            if mid_point == next_midpoint:
                return -1
            
            return binarySearch(spikes, interval, mid_point + 1, end)

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
    while i > -1 or i < len(spikes):

        if inside(spikes[i], interval) and DOWN:
            spike_set.add(i)
            i = i - 1
            
        elif not inside(spikes[i], interval) and DOWN:
            i    = known + 1
            UP   = True
            DOWN = False
            
        elif inside(spikes[i], interval) and UP:
            spike_set.add(i)
            i = i + 1
            
        elif not inside(spikes[i], interval) and UP:
            break

    # Convert set to list, then return.
    return list(spike_set)
    

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


def midpoint(start, end):
    '''
    Description
    -----------
    'midpoint' will calculate midpoint between two points
 
    Input(s)
    --------
    'start'  : int. beginning
    'end'    : int. end

    Output(s)
    --------
    int. midpoint between 'start' and 'end'
    '''
    return int(start + (end - start)/2)
