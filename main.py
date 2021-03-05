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
import matplotlib.pyplot as plt # Matplotlib plotting library


import numpy

from calc_density import calc_density
from calc_density_old import calc_density_old
from calc_density_adapted import calc_density_adapted

def scan_num_cells():
    # constants
    ncells = 20
    npart = 200000
    time_length = 20
    num_time_steps = 50
    L=4*numpy.pi
    constants = [ncells, npart, time_length, num_time_steps,L] # [ncells,npart,time_length,num_time_steps,L]
    
    # variable
    scan_var = 'ncells'#'ncells'/'npart'/'time_length'/'num_time_steps'/'length'
    num_repeats = 10
    start_var = numpy.log2(10)
    stop_var =  numpy.log2(320)
    N = 6
    is_log = True
    this_log_base = 2.0
    
    scan_range(scan_var, num_repeats, start_var, stop_var, N, is_log, constants,this_log_base)

def scan_num_parts():
    # constants
    ncells = 20
    npart = 200000
    time_length = 20
    num_time_steps = 50
    L=4*numpy.pi
    constants = [ncells, npart, time_length, num_time_steps,L] # [ncells,npart,time_length,num_time_steps,L]
    
    # variable
    scan_var = 'npart'#'ncells'/'npart'/'time_length'/'num_time_steps'/'length'
    num_repeats = 10
    start_var = numpy.log10(1000)
    stop_var =  numpy.log10(1000000)
    N = 6
    is_log = True
    this_log_base = 10.0
    
    scan_range(scan_var, num_repeats, start_var, stop_var, N, is_log, constants,this_log_base)

def scan_time_length():
    # constants
    ncells = 20
    npart = 200000
    time_length = 20
    num_time_steps = 50
    L=4*numpy.pi
    constants = [ncells, npart, time_length, num_time_steps,L] # [ncells,npart,time_length,num_time_steps,L]
    
    # variable
    scan_var = 'time_length'#'ncells'/'npart'/'time_length'/'num_time_steps'/'length'
    num_repeats = 10
    start_var = 5
    stop_var = 200
    N = 20
    is_log = False
    
    scan_range(scan_var, num_repeats, start_var, stop_var, N, is_log, constants)

def scan_numTimeSteps():
    # constants
    ncells = 20
    npart = 200000
    time_length = 20
    num_time_steps = 50
    L=4*numpy.pi
    constants = [ncells, npart, time_length, num_time_steps,L]  # [ncells,npart,time_length,num_time_steps,L]
    
    scan_var = 'num_time_steps'#'ncells'/'npart'/'time_length'/'num_time_steps'/'length'
    num_repeats = 10
    start_var = 5
    stop_var =  250
    N = 20
    is_log = False
    this_log_base = 10.0
    
    scan_range(scan_var, num_repeats, start_var, stop_var, N, is_log, constants,this_log_base)

def scan_L():
    # constants
    ncells = 20
    npart = 200000
    time_length = 20
    num_time_steps = 50
    L=4*numpy.pi
    constants = [ncells, npart, time_length, num_time_steps,L]  # [ncells,npart,time_length,num_time_steps,L]
    
    scan_var = 'length'#'ncells'/'npart'/'time_length'/'num_time_steps'/'length'
    num_repeats = 10
    start_var = 0.5*numpy.pi
    stop_var =  10*numpy.pi
    N = 20
    is_log = False
    this_log_base = 10.0
    
    scan_range(scan_var, num_repeats, start_var, stop_var, N, is_log, constants,this_log_base)
    
if __name__ == "__main__":
    #scan_num_cells()
    scan_num_parts()
    #scan_time_length()
    #scan_numTimeSteps()
    #scan_L()
    
    
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
    """
    # constants
    ncells = 20
    npart = 1000
    time_length = 20
    num_time_steps = 50
    constants = [ncells, npart, time_length, num_time_steps] # [ncells,npart,time_length,num_time_steps]
    
    # variables 
    scan_var = 'ncells'#'ncells'/'npart'/'time_length'/'num_time_steps'
    num_repeats = 10
    start_var = numpy.log2(10)
    stop_var =  numpy.log2(320)
    N = 6
    is_log = True
    this_log_base = 2.0
    
    scan_range(scan_var, num_repeats, start_var, stop_var, N, is_log, constants,this_log_base)
    """
    
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