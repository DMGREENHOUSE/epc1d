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
from numpy import linspace, logspace, zeros, std, mean
import matplotlib.pyplot as plt # Matplotlib plotting library

from math import isnan

def scan_ncells(var_val, constants, this_iteration_num):
    var_val = round(var_val)
    this_run = Data_Run(
                    L=constants[4],
                    ncells=int(var_val),
                    npart=int(constants[1]),
                    time_length = constants[2],
                    num_time_steps = int(constants[3]),
                    iteration_num = int(this_iteration_num)
                    )
    return this_run

def scan_npart(var_val, constants, this_iteration_num):
    var_val = round(var_val)
    this_run = Data_Run(
                    L=constants[4],
                    ncells=int(constants[0]),
                    npart=int(var_val),
                    time_length = constants[2],
                    num_time_steps = int(constants[3]),
                    iteration_num = int(this_iteration_num)
                    )
    return this_run

def scan_time_length(var_val, constants, this_iteration_num):
    var_val = round(var_val)
    this_run = Data_Run(
                    L=constants[4],
                    ncells=int(constants[0]),
                    npart=int(constants[1]),
                    time_length = var_val,
                    num_time_steps = int(constants[3]),
                    iteration_num = int(this_iteration_num)
                    )
    return this_run

def scan_num_time_steps(var_val, constants, this_iteration_num):
    var_val = round(var_val)
    this_run = Data_Run(
                    L=constants[4],
                    ncells=int(constants[0]),
                    npart=int(constants[1]),
                    time_length = constants[2],
                    num_time_steps = int(var_val),
                    iteration_num = int(this_iteration_num)
                    )
    return this_run

def scan_length(var_val, constants, this_iteration_num):
    #var_val = round(var_val)
    this_run = Data_Run(
                    L=var_val,
                    ncells=int(constants[0]),
                    npart=int(constants[1]),
                    time_length = constants[2],
                    num_time_steps = int(constants[3]),
                    iteration_num = int(this_iteration_num)
                    )
    return this_run

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
    elif scan_var == 'length':
        func = scan_length
    else:
        print('ERROR: not recognised scan variable')
    return func

def scan_range(scan_var, num_repeats, start, stop, N, is_log, constants, log_base=10):
    arr = None
    if is_log:
        arr = logspace(start, stop, N, base=log_base)
    else:
        arr = linspace(start, stop, N)
    
    func = find_func(scan_var)
    for var_val in arr:
        for i in range(num_repeats):
            new_data_run = func(var_val, constants, i)
            this_summary = new_data_run.get_summary()
            filename = hf.find_filename(new_data_run)
            filepath = r'Data/'+filename
            outfile = open(filepath, 'wb')
            pickle.dump(this_summary,outfile)
            outfile.close()
    
#########################################################
##         For Analysing
#########################################################
def add_filepath(filename):
    filepath = r'Data/'+filename
    return filepath

def perform_analysis(filepath, is_plot_graphs, method):
    infile = open(filepath,'rb')
    new_summary = pickle.load(infile)
    infile.close()
    return Analysis(new_summary, is_plot_graphs, method)
    

def analyse_scan_ncells(var_val, constants, this_iteration_num, type_name, is_plot_graphs, method):
    var_val = round(var_val)
    filename = hf.generate_filename(type_name,
                                    str(constants[4]),
                                    str(int(var_val)),
                                    str(int(constants[1])),
                                    str(constants[2]),
                                    str(int(constants[3])),
                                    str(this_iteration_num))
    return perform_analysis(add_filepath(filename), is_plot_graphs, method)

def analyse_scan_npart(var_val, constants, this_iteration_num, type_name, is_plot_graphs, method):
    var_val = round(var_val)
    filename = hf.generate_filename(type_name,
                                    str(constants[4]),
                                    str(constants[0]),
                                    str(int(var_val)),
                                    str(constants[2]),
                                    str(constants[3]),
                                    str(this_iteration_num))
    return perform_analysis(add_filepath(filename), is_plot_graphs, method)

def analyse_scan_time_length(var_val, constants, this_iteration_num, type_name, is_plot_graphs, method):
    var_val = round(var_val)
    filename = hf.generate_filename(type_name,
                                    str(constants[4]),
                                    str(constants[0]),
                                    str(constants[1]),
                                    str(var_val),
                                    str(constants[3]),
                                    str(this_iteration_num))
    return perform_analysis(add_filepath(filename), is_plot_graphs, method)

def analyse_scan_num_time_steps(var_val, constants, this_iteration_num, type_name, is_plot_graphs, method):
    var_val = round(var_val)
    filename = hf.generate_filename(type_name,
                                    str(constants[4]),
                                    str(constants[0]),
                                    str(constants[1]),
                                    str(constants[2]),
                                    str(int(var_val)),
                                    str(this_iteration_num))
    return perform_analysis(add_filepath(filename), is_plot_graphs, method)

def analyse_scan_length(var_val, constants, this_iteration_num, type_name, is_plot_graphs, method):
    var_val = round(var_val,3)
    filename = hf.generate_filename(type_name,
                                    str(var_val),
                                    str(constants[0]),
                                    str(constants[1]),
                                    str(constants[2]),
                                    str(constants[3]),
                                    str(this_iteration_num))
    return perform_analysis(add_filepath(filename), is_plot_graphs, method)

def find_analyse_func(scan_var):
    func = None
    if scan_var == 'ncells':
        func = analyse_scan_ncells
    elif scan_var == 'npart':
        func = analyse_scan_npart
    elif scan_var == 'time_length':
        func = analyse_scan_time_length
    elif scan_var == 'num_time_steps':
        func = analyse_scan_num_time_steps
    elif scan_var == 'length':
        func = analyse_scan_length
    else:
        print('ERROR: not recognised scan variable')
    return func

def analyse_scan_range(scan_var, num_repeats, start, stop, N, is_log,
                       constants, this_type, is_plot_graphs,
                       is_plot_summary_graph, log_base=10.0):
    arr = None
    if is_log:
        arr = logspace(start, stop, N, base=log_base)
    else:
        arr = linspace(start, stop, N)
    func = find_analyse_func(scan_var)
    scan_range = arr
    
    count = 0
    
    results = zeros((10,N))
    # Results form:
    #0, noise_method
    #1, noise_method_err
    #2, DR_method
    #3, DR_method_err
    #4, time_method
    #5, time_method_err
    #6, freq_method
    #7, freq_method_err
    #8, freq_err_method # the frequency error individually found
    #9, freq_err_method_err
    
    for var_val in arr:
        this_results = [[], # this_noise_method
                        [], #this_DR_method
                        [], #this_time_method, should not change between methods
                        [], #this_freq_method, should not change between methods
                        [] # this_freq_err_method, should not change between methods
                        ]
        for i in range(num_repeats):
            # Method
            method = 'Method 2' # 'Method 1'/'Method 2'/'Method 3'
            new_analysis = func(var_val, constants, i, this_type, is_plot_graphs, method)
            i_results =  new_analysis.get_results()
            for i in range(len(i_results)):
                if i_results[i] is not None and not isnan(i_results[i]):
                    this_results[i].append(i_results[i])
        NULL_VAL = float('nan')
        for i in range(len(this_results)):
            if len(this_results[i]) == 0:
                results[i][count] = NULL_VAL
            else:
                index = int(2*i)
                results[index][count] = mean(this_results[i])
                results[index+1][count] = std(this_results[i])
        count+=1
    if is_plot_summary_graph:
        print(scan_range)
        print(results[0])
        print(results[1])
        plt.errorbar(scan_range, results[0], yerr=results[1], marker='x', label = 'log method')
        plt.ioff() # This so that the windows stay open
        plt.legend()
        plt.show()
    return [scan_range, results[0], results[1], results[4], results[5], results[6], results[7], results[2], results[3]]
    