#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 18:27:09 2021

@author: dmg530
"""
from Data_Run import Data_Run
from Analysis import Analysis
import pickle
import helper_functions as hf
import time


if __name__ == "__main__":
    type_name = 'LD' # 'LD' or '2SI'
    length = '12.566' # length as a string
    num_cells = '100' # number of cells as a string
    num_parts = '10000' # number of particles as a string
    time_length = '50' # length of time as a string
    time_steps = '100' # number of timesteps as a string
    iteration_num = '1' # iteration number as a string
    
    filename = hf.generate_filename(type_name, length, num_cells, num_parts,
                                    time_length, time_steps, iteration_num)
    filepath = r'Data/'+filename
    
    # input the object
    infile = open(filepath,'rb')
    new_summary = pickle.load(infile)
    infile.close()
    new_analysis = Analysis(new_summary)