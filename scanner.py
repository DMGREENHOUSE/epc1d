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

def perform_analysis(filepath, is_plot_graphs):
    infile = open(filepath,'rb')
    new_summary = pickle.load(infile)
    infile.close()
    return Analysis(new_summary, is_plot_graphs)
    

def analyse_scan_ncells(var_val, constants, this_iteration_num, type_name, is_plot_graphs):
    var_val = round(var_val)
    filename = hf.generate_filename(type_name,
                                    str(constants[4]),
                                    str(int(var_val)),
                                    str(int(constants[1])),
                                    str(constants[2]),
                                    str(int(constants[3])),
                                    str(this_iteration_num))
    return perform_analysis(add_filepath(filename), is_plot_graphs)

def analyse_scan_npart(var_val, constants, this_iteration_num, type_name, is_plot_graphs):
    var_val = round(var_val)
    filename = hf.generate_filename(type_name,
                                    str(constants[4]),
                                    str(constants[0]),
                                    str(int(var_val)),
                                    str(constants[2]),
                                    str(constants[3]),
                                    str(this_iteration_num))
    return perform_analysis(add_filepath(filename), is_plot_graphs)

def analyse_scan_time_length(var_val, constants, this_iteration_num, type_name, is_plot_graphs):
    filename = hf.generate_filename(type_name,
                                    str(constants[4]),
                                    str(constants[0]),
                                    str(constants[1]),
                                    str(var_val),
                                    str(constants[3]),
                                    str(this_iteration_num))
    return perform_analysis(add_filepath(filename), is_plot_graphs)

def analyse_scan_num_time_steps(var_val, constants, this_iteration_num, type_name, is_plot_graphs):
    var_val = round(var_val)
    filename = hf.generate_filename(type_name,
                                    str(constants[4]),
                                    str(constants[0]),
                                    str(constants[1]),
                                    str(constants[2]),
                                    str(int(var_val)),
                                    str(this_iteration_num))
    return perform_analysis(add_filepath(filename), is_plot_graphs)

def analyse_scan_length(var_val, constants, this_iteration_num, type_name, is_plot_graphs):
    var_val = round(var_val,3)
    filename = hf.generate_filename(type_name,
                                    str(var_val),
                                    str(constants[0]),
                                    str(constants[1]),
                                    str(constants[2]),
                                    str(constants[3]),
                                    str(this_iteration_num))
    return perform_analysis(add_filepath(filename), is_plot_graphs)

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
    
    noise_log_method = zeros([N])
    noise_log_method_err = zeros([N])
    DR_log_method = zeros([N])
    DR_log_method_err = zeros([N])
    
    noise_nolog_method = zeros([N])
    noise_nolog_method_err = zeros([N])
    DR_nolog_method = zeros([N])
    DR_nolog_method_err = zeros([N])
    
    noise_noabs_method = zeros([N])
    noise_noabs_method_err = zeros([N])
    DR_noabs_method = zeros([N])
    DR_noabs_method_err = zeros([N])
    count = 0
    times=[None]*N
    time_errs=[None]*N
    for var_val in arr:
        this_noise_log_method = []
        this_DR_log_method = []
        this_noise_nolog_method = []
        this_DR_nolog_method = []
        this_noise_noabs_method = []
        this_DR_noabs_method = []
        
        this_time=[]
        for i in range(num_repeats):
            new_analysis = func(var_val, constants, i, this_type, is_plot_graphs)
            results_log =  new_analysis.get_results_log()
            results_no_log =  new_analysis.get_results_nolog()
            results_no_abs =  new_analysis.get_results_noabs()
            results_time =  new_analysis.get_results_time()
            results_frequency =  new_analysis.get_results_frequency()
            results_frequency_err =  new_analysis.get_results_frequency_err()
            print('LOOK: ')
            print(results_frequency, results_frequency_err)
            if results_log[0] is not None:
                this_noise_log_method.append(results_log[0])
            if results_log[1] is not None:
                this_DR_log_method.append(results_log[1])
                
            if results_no_log[0] is not None:
                this_noise_nolog_method.append(results_no_log[0])
            if results_no_log[1] is not None:
                this_DR_nolog_method.append(results_no_log[1])
                
            if results_no_abs[0] is not None:
                this_noise_noabs_method.append(results_no_abs[0])
            if results_no_abs[1] is not None:
                this_DR_noabs_method.append(results_no_abs[1])
                
            if results_time is not None:
                this_time.append(results_time)
        NULL_VAL = -1
        if len(this_noise_log_method) == 0:
            noise_log_method[count] = NULL_VAL
        else:
            noise_log_method[count] = mean(this_noise_log_method)
            noise_log_method_err[count] = std(this_noise_log_method)
        
        if len(this_DR_log_method) == 0:
            DR_log_method[count] = NULL_VAL
        else:
            DR_log_method[count] = mean(this_DR_log_method)
            DR_log_method_err[count] = std(this_DR_log_method)
        
        if len(this_noise_nolog_method) == 0:
            noise_nolog_method[count] = NULL_VAL
        else:
            noise_nolog_method[count] = mean(this_noise_nolog_method)
            noise_nolog_method_err[count] = std(this_noise_nolog_method)
        
        if len(this_DR_nolog_method) == 0:
            DR_nolog_method[count] = NULL_VAL
        else:
            DR_nolog_method[count] = mean(this_DR_nolog_method)
            DR_nolog_method_err[count] = std(this_DR_nolog_method)
        
        if len(this_noise_noabs_method) == 0:
            noise_noabs_method[count] = NULL_VAL
        else:
            noise_noabs_method[count] = mean(this_noise_noabs_method)
            noise_noabs_method_err[count] = std(this_noise_noabs_method)
        
        if len(this_DR_noabs_method) == 0:
            DR_noabs_method[count] = NULL_VAL
        else:
            DR_noabs_method[count] = mean(this_DR_noabs_method)
            DR_noabs_method_err[count] = std(this_DR_noabs_method)
        
        if len(this_time) == 0:
            times[count] = NULL_VAL
            time_errs[count] = NULL_VAL
        else:
            times[count] = mean(this_time)
            time_errs[count] = std(this_time)
        count+=1
    if is_plot_summary_graph:
        print(scan_range)
        print(noise_nolog_method)
        print(noise_nolog_method_err)
        plt.errorbar(scan_range, noise_log_method, yerr=noise_log_method_err, marker='x', label = 'log method')
        plt.errorbar(scan_range, noise_nolog_method, yerr=noise_nolog_method_err, marker='x', label = 'no log method')
        plt.errorbar(scan_range, noise_noabs_method, yerr=noise_noabs_method_err, marker='x', label = 'no log, no abs method')
        plt.ioff() # This so that the windows stay open
        plt.legend()
        plt.show()
    return [scan_range, noise_nolog_method, noise_nolog_method_err, times, time_errs]