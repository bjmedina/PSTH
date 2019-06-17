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
    table['NaN'] = np.zeros((num_bins, 2))
    table['NaN'][:, 0] = bins
    
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
