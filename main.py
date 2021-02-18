# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 13:49:06 2021

@author: dmg530
"""

from Data_Run import Data_Run
from Analysis import Analysis
import pickle
import helper_functions as hf
import time
from scanner import scan_range

import numpy

from calc_density import calc_density
from calc_density_old import calc_density_old
from calc_density_adapted import calc_density_adapted

if __name__ == "__main__":
    """
    N = 500000
    L = 10
    pos = numpy.random.uniform(low=0, high=L, size=(N,))
    ncells = 20
    print(pos)
    start = time.time()
    
    calc_density_old_results = calc_density_old(pos, ncells, L)
    calc_density_old_time = time.time()
    
    calc_density_results = calc_density(pos, ncells, L)
    calc_density_time = time.time()
    
    calc_density_adapted_results = calc_density_adapted(pos, ncells, L)
    calc_density_adapted_time = time.time()
        
    
    print('calc_density_old ', calc_density_old_time - start)
    print('calc_density ', calc_density_time - calc_density_old_time)
    print('calc_density_adapted ', calc_density_adapted_time - calc_density_time)
    
    print(calc_density_results)
    print(calc_density_old_results)
    print(calc_density_adapted_results)
    
    """
    scan_var = 'npart'#'ncells'/'npart'/'time_length'/'num_time_steps'
    num_repeats = 1
    start_var = 1000
    stop_var = 1000
    N = 1
    is_log = False
    constants = [100, 100, 50, 100] # [ncells,npart,time_length,num_time_steps]
    
    scan_range(scan_var, num_repeats, start_var, stop_var, N, is_log, constants)
    """
    start = time.time()
    new_data_run = Data_Run(ncells=10, npart=1000,time_length = 50, num_time_steps = 100)
    end = time.time()
    time_taken = end-start
    print('TIME TAKEN: ', time_taken)
    this_summary = new_data_run.get_summary()
    #new_analysis = Analysis(this_summary)
    filename = hf.find_filename(new_data_run)
    filepath = r'Data/'+filename
    
    
    outfile = open(filepath, 'wb')
    pickle.dump(this_summary,outfile)
    outfile.close()
    
    constants = [] # including scan var
    scan_var = 1# variable to scan over
    scan_range(scan_var, num_repeats, 1, 10, 10, False, constants)
    scan_range(scan_var, num_repeats, 1, 4, 4, True, constants)
    """"""
    
    # input the object
    infile = open(filename,'rb')
    new_summary = pickle.load(infile)
    infile.close()
    new_analysis = Analysis(new_summary)
    """