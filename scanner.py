#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 20:09:24 2021

@author: dmg530
"""

from Data_Run import Data_Run
from Analysis import Analysis
import pickle
import helper_functions as hf
import time
from numpy import linspace, logspace

def scan_ncells(var_val, constants, this_iteration_num):
    print(var_val, constants[1], constants[2], constants[3], this_iteration_num)
    return Data_Run(ncells=int(var_val),
                    npart=int(constants[1]),
                    time_length = constants[2],
                    num_time_steps = int(constants[3]),
                    iteration_num = int(this_iteration_num)
                    )

def scan_npart(var_val, constants, this_iteration_num):
    return Data_Run(ncells=int(constants[0]),
                    npart=int(var_val),
                    time_length = constants[2],
                    num_time_steps = int(constants[3]),
                    iteration_num = int(this_iteration_num)
                    )

def scan_time_length(var_val, constants, this_iteration_num):
    return Data_Run(ncells=int(constants[0]),
                    npart=int(constants[1]),
                    time_length = var_val,
                    num_time_steps = int(constants[3]),
                    iteration_num = int(this_iteration_num)
                    )

def scan_num_time_steps(var_val, constants, this_iteration_num):
    return Data_Run(ncells=int(constants[0]),
                    npart=int(constants[1]),
                    time_length = constants[2],
                    num_time_steps = int(var_val),
                    iteration_num = int(this_iteration_num)
                    )

def find_func(scan_var):
    func = None
    if scan_var == 'ncells':
        func = scan_ncells
    elif scan_var == 'npart':
        func = scan_npart
    elif scan_var == 'time_length':
        func = scan_time_length
    elif scan_var == 'num_time_steps':
        func = scan_num_time_steps
    else:
        print('ERROR: not recognised scan variable')
    return func

def scan_range(scan_var, num_repeats, start, stop, N, is_log, constants):
    arr = None
    if is_log:
        arr = logspace(start, stop, N, base=10.0)
    else:
        arr = linspace(start, stop, N)
    print(arr)
    func = find_func(scan_var)
    for i in range(num_repeats):
        for var_val in arr:
            new_data_run = func(var_val, constants, i)
            this_summary = new_data_run.get_summary()
            filename = hf.find_filename(new_data_run)
            filepath = r'Data/'+filename
            outfile = open(filepath, 'wb')
            pickle.dump(this_summary,outfile)
            outfile.close()
    