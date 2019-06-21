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
import pickle
########################

width = 1

### This is for the 'drifting gratings' stimulus
# All possible temporal frequencies for the stimulus
temp_freqs   = [1, 2, 4, 8, 15]

# All possible orientations of stimulus (angles)
orientations = [i*45 for i in range(8)]
    
start   = 0
end     = 2000
mapping = {'probeA': 'AM',
        'probeB': 'PM',
        'probeC': 'V1',
        'probeD': 'LM',
        'probeE': 'AL',
        'probeF': 'RL'}

###

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
        # I think this is how it works... But it's the idea
        # Has to be a little more complicated than this
        # You need to put it in the right bin

        bn = insertToBin(spike, end)
        #print("config: " + (config) + "\t spike: " + str(spike) + "\t bin: " + str(bn))
        self.__table[config][bn] += 1
        #if(self.__table[config][bn] > 1):
            #print("BIN: " + str(self.__table[config][bn]))

def makeTable():
    '''
    Description
    -----------
    'makeTable' creates a dictionary to keep track of time bins for each possible orientation of stimulus. One for each cell.
    
    Output(s)
    ---------
    'table': dict. Dictionary that contains orientation combination as key and all cells that are in V.
    '''

    bins = np.linspace(start, end, int( (end - start)/width + 1 )) 
    
    # In this table, each key is a different configuration of the stimulus
    # and each row corresponds to spikes in a time bin.
    table        = {}

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
    'start'  : int. beginning
    'end'    : int. end


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

                           
def insertToBin(spiketime, end):
    
    ## TODO
    # I think best thing to do is to wrong 3 down to next multiple of timebins
    # PROBABLY have to subtract "start"
    # This doesn't support spikes that are not in the range.....
    
    idx = int( (spiketime - (spiketime % width)) / width )
    
    if( idx > end ): 
        #print("spiketime " + str(spiketime) + "\tidx " + str(idx))
        idx = end

    return idx 



def saveData(MOUSE_ID):
    # open NWB files using h5py library
    
    ###### UPDATE PATH #################################
    DIRECTORY = '/Users/bjm/Documents/CMU/Research/data'
        
    # Changes depending on the trial.
    start     = 0 #in second
    end       = 2000 #in seconds
    
    # Get file from directory
    spikes_nwb_file = os.path.join(DIRECTORY, 'mouse' + MOUSE_ID + '.spikes.nwb')
    nwb = h5.File(spikes_nwb_file, 'r')
    
    # names of probes
    names = nwb['processing'].keys()
    
    # time stamps ( this never changes )
    # This is SPECIFICALLY for the 'drifting_gratings_2' stimulus
    timestamps  = nwb['stimulus']['presentation']['drifting_gratings_2']['timestamps'].value
    stim_orient = nwb['stimulus']['presentation']['drifting_gratings_2']['data'].value
    
    
    ## Adding spikes
    for probe_name in names:
        # Get all cells that are in V for every probe
        #print(probe_name)
        probe = Probe(nwb, probe_name)
        
        # Going to want to save this information later.
        filename = MOUSE_ID + "_" + probe_name
        
        # ...get every cell. Then...
        cells = probe.getCellList()
        
        # ... for every cell...
        for cell in cells:
            
            #print("\tNeuron:\t" + str(cell))
    
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
                            curr_cell.addSpike(config, stim_spike, end )
                            
        print("Saving to " + filename)
    
        with open(filename, 'wb') as f:
                pickle.dump(probe, f)
                
