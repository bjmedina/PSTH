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
        Method returns dictionary of cells
        
        Output(s)
        ---------
        Dictionary of cell 'cell_number'
        '''

        return self.__cells.keys()


class Cell:
    
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

    def getSpikes(self, orientation):
        '''
        Description
        -----------
        Method returns table for given cell
        
        Input(s)
        --------
        'orientation': string. key of dictionary.
        
        Output(s)
        ---------
        table at certain orientation
        '''
        return self.__table[orientation]

    def addSpike(self, orientation, spike):
        # I think this is how it works... But it's the idea 
        self.__table[orientation].append(spike)

        

    ### Add function for filling table
    
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


def makeTable():
    '''
    Description
    -----------
    'makeTable' creates a dictionary to keep track of time bins for each possible orientation of stimulus. One for each cell.
    
    Output(s)
    ---------
    'table': dict. Dictionary that contains orientation combination as key and all cells that are in V.
    '''
    
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
            table[config] = []
            
    # Empty images
    table['nan_nan'] = []
     
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
            
            return binarySearch(spikes, interval, start, mid_point-1)

        elif spikes[mid_point] < interval[0]:
            
            next_midpoint = midpoint(mid_point+1, end)
            
            # If this is true, then we're going to hit a recursion error....
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
    while i > -1 or i < len(spikes):

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


                           
