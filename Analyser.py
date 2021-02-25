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
from scanner import analyse_scan_range
import numpy
import matplotlib.pyplot as plt # Matplotlib plotting library
from scipy import interpolate, optimize


COLOURS=['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
        'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

def lin_func2(fit_vals,t):
    return fit_vals[0]*(t)+fit_vals[1]
        
def pow_func2(fit_vals,t):
    return fit_vals[0]*pow((t),fit_vals[1])

def perform_fit2(xs, ys,variables = [1,1], func=pow_func2):
    fit_vals, cov_vals = optimize.curve_fit(
                    lambda t,a,b: func([a,b],t),
                    xs,  ys,  p0=(variables[0], variables[1]),
                                        maxfev=1000)
    errs = [cov_vals[0][0], cov_vals[1][1]]
    return fit_vals, errs
        
def pow_func3(fit_vals,t):
    return fit_vals[0]*pow((t),fit_vals[1])+fit_vals[2]
        
def quad_func3(fit_vals,t):
    return fit_vals[0]*pow((t-fit_vals[1]),2)+fit_vals[2]

def perform_fit3(xs, ys,variables = [1,1,0], func=pow_func3):
    fit_vals, cov_vals = optimize.curve_fit(
                    lambda t,a,b,c: func([a,b,c],t),
                    xs,  ys,  p0=(variables[0], variables[1],variables[2]),
                                        maxfev=1000)
    errs = [cov_vals[0][0], cov_vals[1][1], cov_vals[2][2]]
    return fit_vals, errs

def perform_single_analysis():
    type_name = 'LD' # 'LD' or '2SI'
    length = 12.566 # length as a string
    num_cells = '100' # number of cells as a string
    num_parts = '100000' # number of particles as a string
    time_length = '50' # length of time as a string
    time_steps = '100' # number of timesteps as a string
    iteration_num = '0' # iteration number as a string
    
    filename = hf.generate_filename(type_name, length, num_cells, num_parts,
                                    time_length, time_steps, iteration_num)
    filepath = r'Data/'+filename
    
    # input the object
    infile = open(filepath,'rb')
    new_summary = pickle.load(infile)
    infile.close()
    new_analysis = Analysis(new_summary)
    
def perform_single_scan(scan_var, num_repeats, start_var, stop_var, N,
                       is_log, constants, type_name,
                       is_plot_graphs, is_plot_summary_graph, this_log_base=10):
    return analyse_scan_range(scan_var, num_repeats, start_var, stop_var, N,
                       is_log, constants, type_name,
                       is_plot_graphs, is_plot_summary_graph, log_base=this_log_base)

def plot_scan_noise_vs_numCells():
    """ Find and plot scan for noise floor across the number of cells
    Performed for various numbers of particles
    """
    
    # constants
    ncells = 10
    npart = 100
    time_length = 50
    num_time_steps = 100
    length = 12.566 # length as a string
    constants = [ncells, npart, time_length, num_time_steps,length]
    
    # type
    type_name = 'LD' # 'LD' or '2SI'
    
    # variable
    scan_var = 'ncells'#'ncells'/'npart'/'time_length'/'num_time_steps'
    num_repeats = 10
    start = 10
    stop = 320
    start_var = numpy.log2(start)
    stop_var =  numpy.log2(stop)
    N = 6
    is_log = True
    this_log_base = 2.0
    
    # general 
    is_plot_graphs = False
    is_plot_summary_graph = False
    
    # secondary scan
    secondary_vals = [10,100,1000,10000,100000]
    num_sec_vals = len(secondary_vals)
    results = [None]*num_sec_vals
    
    METHOD_NAME = 'Method 2' # for legend
    SECONDARY_SCAN_NAME = 'Particles' # for legend
    
    for i in range(num_sec_vals):
        constants[1] = secondary_vals[i]
        results[i] = analyse_scan_range(scan_var, num_repeats, start_var, stop_var, N,
                       is_log, constants, type_name,
                       is_plot_graphs, is_plot_summary_graph, log_base=this_log_base)
        this_label = METHOD_NAME+', {} '.format(secondary_vals[i])+SECONDARY_SCAN_NAME
        plt.errorbar(results[i][0], results[i][1], yerr=results[i][2], color = COLOURS[i], marker='x', label = this_label)
    
        # Fit
        fit_vals, fit_square_errs = perform_fit2(results[i][0], results[i][1], variables=[1,0,0.1], func=lin_func2)
        num_fit_points = 100
        fit_xs = numpy.linspace(start, stop, num_fit_points)
        fit_ys = lin_func2(fit_vals, fit_xs)
        fit_label = r'%.0f Particles Linear Fit: $(%.2g\pm %.1g)(x)+(%.2g\pm %.1g)$' % (
                secondary_vals[i],
                fit_vals[0],
                numpy.sqrt(fit_square_errs[0]),
                fit_vals[1],
                numpy.sqrt(fit_square_errs[1])
                )
        
        plt.plot(fit_xs, fit_ys, linestyle='--', color = COLOURS[i], label=fit_label)
    
    plt.ioff() # This so that the windows stay open
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('Noise Floor Value')
    plt.xlabel('Number of Cells')
    plt.show()

def plot_scan_noise_vs_numPart():
    """ Find and plot scan for noise floor across the number of particles
    Performed for various numbers of cells
    """
    # constants
    ncells = 10
    npart = 100
    time_length = 50
    num_time_steps = 100
    length = 12.566 # length as a string
    constants = [ncells, npart, time_length, num_time_steps,length]
    
    # type
    type_name = 'LD' # 'LD' or '2SI'
    
    # variable
    scan_var = 'npart'#'ncells'/'npart'/'time_length'/'num_time_steps'
    num_repeats = 10
    start = 10
    stop = 100000
    start_var = numpy.log10(start)
    stop_var =  numpy.log10(stop)
    
    
    N = 5
    is_log = True
    
    # general 
    is_plot_graphs = False
    is_plot_summary_graph = False
    
    # secondary scan
    secondary_vals = [10,20,40,80,160,320]
    num_sec_vals = len(secondary_vals)
    results = [None]*num_sec_vals
    
    METHOD_NAME = 'Method 2' # for legend
    SECONDARY_SCAN_NAME = 'Cells' # for legend
    
    for i in range(num_sec_vals):
        constants[0] = secondary_vals[i]
        results[i] = analyse_scan_range(scan_var, num_repeats, start_var, stop_var, N,
                       is_log, constants, type_name,
                       is_plot_graphs, is_plot_summary_graph)
        this_label = METHOD_NAME+', {} '.format(secondary_vals[i])+SECONDARY_SCAN_NAME
        print(results[i][1])
        plt.errorbar(results[i][0], results[i][1], yerr=results[i][2], color = COLOURS[i], marker='x', label = this_label)
    
        # Fit
        func = pow_func3
        fit_vals, fit_square_errs = perform_fit3(results[i][0], results[i][1], variables=[1,-0.5,1])
        num_fit_points = 100
        fit_xs = numpy.linspace(start, stop, num_fit_points)
        fit_ys = func(fit_vals, fit_xs)
        fit_label = r'%.0f Cells Linear Fit: $(%.2g\pm %.1g)(x)^{(%.2g\pm %.1g)}+(%.2g\pm %.1g)$' % (
                secondary_vals[i],
                fit_vals[0],
                numpy.sqrt(fit_square_errs[0]),
                fit_vals[1],
                numpy.sqrt(fit_square_errs[1]),
                fit_vals[2],
                numpy.sqrt(fit_square_errs[2])
                )
        
        plt.plot(fit_xs, fit_ys, linestyle='--',color = COLOURS[i], label=fit_label)
    plt.ioff() # This so that the windows stay open
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    #plt.xscale('linear')
    #plt.yscale('linear')
    plt.ylabel('Noise Floor Value')
    plt.xlabel('Number of Particles')
    plt.show()

def plot_scan_noise_vs_timeLength():
    """ Find and plot scan for noise floor across the number of particles
    Performed for var
    """
    # constants
    ncells = 20
    npart = 1000
    time_length = 20
    num_time_steps = 50
    length = 12.566 # length as a string
    constants = [ncells, npart, time_length, num_time_steps,length]
    
    # type
    type_name = 'LD' # 'LD' or '2SI'
    
    # variable
    scan_var = 'time_length'#'ncells'/'npart'/'time_length'/'num_time_steps'
    num_repeats = 10
    start_var = 5
    stop_var = 200
    N = 20
    is_log = False
    
    # general 
    is_plot_graphs = False
    is_plot_summary_graph = False
    
    results = analyse_scan_range(scan_var, num_repeats, start_var, stop_var, N,
                       is_log, constants, type_name,
                       is_plot_graphs, is_plot_summary_graph)
    
    # Fit
    this_func=lin_func2
    fit_vals, fit_square_errs = perform_fit2(results[0][1:], results[1][1:],func=this_func)
    num_fit_points = 100
    fit_xs = numpy.linspace(15, stop_var, num_fit_points)
    fit_ys = this_func(fit_vals, fit_xs)
    fit_label = r'Linear Fit: $(%.2g\pm %.1g)(x)+(%.2g\pm %.1g)$' % (
            fit_vals[0],
            numpy.sqrt(fit_square_errs[0]),
            fit_vals[1],
            numpy.sqrt(fit_square_errs[1]))
    plt.plot(fit_xs, fit_ys, label=fit_label)
    
    
    this_label = hf.generate_filename(
            str(type_name), str(length), str(ncells), str(npart),
            'VAR', str(num_time_steps), 'AVG.')
    plt.errorbar(results[0], results[1], yerr=results[2], color = COLOURS[0], marker='x', label = this_label)
    plt.ioff() # This so that the windows stay open
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    #plt.xscale('linear')
    #plt.yscale('linear')
    plt.ylabel('Noise Floor Value')
    plt.xlabel('Simulation Output Time Length')
    plt.show()

def plot_scan_noise_vs_numTimeSteps():
    """ Find and plot scan for noise floor across the number of particles
    Performed for var
    """
    # constants
    ncells = 20
    npart = 1000
    time_length = 20
    num_time_steps = 50
    length = 12.566 # length as a string
    constants = [ncells, npart, time_length, num_time_steps, length]
    
    # type
    type_name = 'LD' # 'LD' or '2SI'
    
    # variable
    scan_var = 'num_time_steps'#'ncells'/'npart'/'time_length'/'num_time_steps'/'length'
    num_repeats = 10
    start_var = 5
    stop_var =  250
    N = 20
    is_log = False
    
    # general 
    is_plot_graphs = False
    is_plot_summary_graph = False
    
    results = analyse_scan_range(scan_var, num_repeats, start_var, stop_var, N,
                       is_log, constants, type_name,
                       is_plot_graphs, is_plot_summary_graph)
    
    # Fit
    this_func=lin_func2
    fit_vals, fit_square_errs = perform_fit2(results[0][1:], results[1][1:],func=this_func)
    num_fit_points = 100
    fit_xs = numpy.linspace(15, stop_var, num_fit_points)
    fit_ys = this_func(fit_vals, fit_xs)
    fit_label = r'Linear Fit: $(%.2g\pm %.1g)(x)+(%.2g\pm %.1g)$' % (
            fit_vals[0],
            numpy.sqrt(fit_square_errs[0]),
            fit_vals[1],
            numpy.sqrt(fit_square_errs[1]))
    plt.plot(fit_xs, fit_ys, label=fit_label)
    
    
    this_label = hf.generate_filename(
            str(type_name), str(length), 'VAR', str(npart),
            str(time_length), str(num_time_steps), 'AVG:{}'.format(num_repeats))
    plt.errorbar(results[0], results[1], yerr=results[2], color = COLOURS[0], marker='x', label = this_label)
    plt.ioff() # This so that the windows stay open
    plt.legend()
    plt.xscale('linear')
    plt.yscale('linear')
    plt.ylabel('Noise Floor Value')
    plt.xlabel('Number of Timesteps')
    plt.show()

def plot_scan_noise_vs_length():
    """ Find and plot scan for noise floor across the number of particles
    Performed for var
    """
    # constants
    ncells = 20
    npart = 1000
    time_length = 20
    num_time_steps = 50
    length = 12.566 # length as a string
    constants = [ncells, npart, time_length, num_time_steps, length]
    
    # type
    type_name = 'LD' # 'LD' or '2SI'
    
    # variable
    scan_var = 'length'#'ncells'/'npart'/'time_length'/'num_time_steps'/'length'
    num_repeats = 10
    start_var = 0.5*numpy.pi
    stop_var =  10*numpy.pi
    N = 20
    is_log = False
    
    # general 
    is_plot_graphs = False
    is_plot_summary_graph = False
    
    results = analyse_scan_range(scan_var, num_repeats, start_var, stop_var, N,
                       is_log, constants, type_name,
                       is_plot_graphs, is_plot_summary_graph)
    
    # Fit
    this_func=quad_func3
    fit_vals, fit_square_errs = perform_fit3(results[0], results[1],func=this_func,variables=[0.1,4,0.05])
    num_fit_points = 100
    fit_xs = numpy.linspace(start_var, stop_var, num_fit_points)
    fit_ys = this_func(fit_vals, fit_xs)
    fit_label = r'Quadratic Fit: $(%.2g\pm %.1g)(x-(%.2g\pm %.1g))^{2}+(%.2g\pm %.1g)$' % (
            fit_vals[0],
            numpy.sqrt(fit_square_errs[0]),
            fit_vals[1],
            numpy.sqrt(fit_square_errs[1]),
            fit_vals[2],
            numpy.sqrt(fit_square_errs[2]))
    plt.plot(fit_xs, fit_ys, linestyle='--', label=fit_label)
    
    
    this_label = hf.generate_filename(
            str(type_name), str(length), 'VAR', str(npart),
            str(time_length), str(num_time_steps), 'AVG:{}'.format(num_repeats))
    plt.errorbar(results[0], results[1], yerr=results[2],color = COLOURS[0], marker='x', label = this_label)
    plt.axvline(x=4*numpy.pi, linestyle='dotted',label='$4\pi$')
    plt.ioff() # This so that the windows stay open
    plt.legend()
    plt.xscale('linear')
    plt.yscale('linear')
    plt.ylabel('Noise Floor Value')
    plt.xlabel('Length')
    plt.show()

#############################################################
##### TIME PLOTS
#############################################################
def plot_time_vs_timeLength():
    # constants
    ncells = 20
    npart = 1000
    time_length = 20
    num_time_steps = 50
    length = 12.566 # length as a string
    constants = [ncells, npart, time_length, num_time_steps,length]
    
    # type
    type_name = 'LD' # 'LD' or '2SI'
    
    # variable
    scan_var = 'time_length'#'ncells'/'npart'/'time_length'/'num_time_steps'
    num_repeats = 10
    start_var = 5
    stop_var = 200
    N = 20
    is_log = False
    
    # general 
    is_plot_graphs = False
    is_plot_summary_graph = False
    results = analyse_scan_range(scan_var, num_repeats, start_var, stop_var, N,
                       is_log, constants, type_name,
                       is_plot_graphs, is_plot_summary_graph)
    
    # Fit
    fit_vals, fit_square_errs = perform_fit2(results[0], results[3])
    num_fit_points = 100
    fit_xs = numpy.linspace(start_var, stop_var, num_fit_points)
    fit_ys = pow_func2(fit_vals, fit_xs)
    fit_label = r'Power Fit: $(%.2g\pm %.1g)(x)^{%.2g\pm %.1g}$' % (
            fit_vals[0],
            numpy.sqrt(fit_square_errs[0]),
            fit_vals[1],
            numpy.sqrt(fit_square_errs[1]))
    plt.plot(fit_xs, fit_ys, label=fit_label)
    plt.yscale('linear')
    
    this_label = hf.generate_filename(
            str(type_name), str(length), str(ncells), str(npart),
            'VAR', str(num_time_steps), 'AVG.')
    plt.errorbar(results[0], results[3], yerr=results[4], marker='x', label = this_label)
    plt.ylabel('Time Taken [s]')
    plt.xlabel('Simulation Output Time Length')
    plt.ioff() # This so that the windows stay open
    plt.legend()
    plt.show()

def plot_time_vs_numCells():
    # constants
    ncells = 20
    npart = 1000
    time_length = 20
    num_time_steps = 50
    length = 12.566 # length as a string
    constants = [ncells, npart, time_length, num_time_steps,length]
    
    # type
    type_name = 'LD' # 'LD' or '2SI'
    
    # variable
    scan_var = 'ncells'#'ncells'/'npart'/'time_length'/'num_time_steps'
    num_repeats = 10
    start_var = numpy.log2(10)
    stop_var =  numpy.log2(320)
    N = 6
    is_log = True
    this_log_base = 2.0
    
    # general 
    is_plot_graphs = False
    is_plot_summary_graph = False
    results = analyse_scan_range(scan_var, num_repeats, start_var, stop_var, N,
                       is_log, constants, type_name,
                       is_plot_graphs, is_plot_summary_graph, this_log_base)
    
    this_label = hf.generate_filename(
            str(type_name), str(length), 'VAR', str(npart),
            str(time_length), str(num_time_steps), 'AVG:{}'.format(num_repeats))
    
    # Fit
    fit_vals, fit_square_errs = perform_fit2(results[0], results[3])
    num_fit_points = 100
    fit_xs = numpy.logspace(start_var, stop_var, num_fit_points, base=this_log_base)
    fit_ys = pow_func2(fit_vals, fit_xs)
    fit_label = r'Power Fit: $(%.2g\pm %.1g)(x)^{%.2g\pm %.1g}$' % (
            fit_vals[0],
            numpy.sqrt(fit_square_errs[0]),
            fit_vals[1],
            numpy.sqrt(fit_square_errs[1]))
    
    plt.yscale('linear')
    plt.plot(fit_xs, fit_ys, label=fit_label)
    
    plt.errorbar(results[0], results[3], yerr=results[4], marker='x', label = this_label)
    plt.ylabel('Time Taken [s]')
    plt.xlabel('Number of Cells')
    plt.ioff() # This so that the windows stay open
    plt.legend()
    plt.show()

def plot_time_vs_numParts():
    # constants
    ncells = 20
    npart = 1000
    time_length = 20
    num_time_steps = 50
    length = 12.566 # length as a string
    constants = [ncells, npart, time_length, num_time_steps, length]
    
    # type
    type_name = 'LD' # 'LD' or '2SI'
    
    # variable
    scan_var = 'npart'#'ncells'/'npart'/'time_length'/'num_time_steps'
    num_repeats = 10
    start_var = numpy.log10(10)
    stop_var =  numpy.log10(1000000)
    N = 6
    is_log = True
    this_log_base = 10.0
    
    # general 
    is_plot_graphs = False
    is_plot_summary_graph = False
    results = analyse_scan_range(scan_var, num_repeats, start_var, stop_var, N,
                       is_log, constants, type_name,
                       is_plot_graphs, is_plot_summary_graph, this_log_base)
    
    # Fit
    fit_vals, fit_square_errs = perform_fit3(results[0], results[3])
    num_fit_points = 100
    fit_xs = numpy.logspace(start_var, stop_var, num_fit_points, base=this_log_base)
    fit_ys = pow_func3(fit_vals, fit_xs)
    fit_label = r'Power Fit: $(%.2g\pm %.1g)(x)^{(%.2g\pm %.1g)}+(%.2g\pm %.1g)$' % (
            fit_vals[0],
            numpy.sqrt(fit_square_errs[0]),
            fit_vals[1],
            numpy.sqrt(fit_square_errs[1]),
            fit_vals[2],
            numpy.sqrt(fit_square_errs[2]))
    
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(fit_xs, fit_ys, label=fit_label)
    
    this_label = hf.generate_filename(
            str(type_name), str(length), 'VAR', str(npart),
            str(time_length), str(num_time_steps), 'AVG:{}'.format(num_repeats))
    plt.errorbar(results[0], results[3], yerr=results[4], marker='x', label = this_label)
    plt.ylabel('Time Taken [s]')
    plt.xlabel('Number of Particles')
    plt.ioff() # This so that the windows stay open
    plt.legend()
    plt.show()

def plot_time_vs_numTimeSteps():
    # constants
    ncells = 20
    npart = 1000
    time_length = 20
    num_time_steps = 50
    length = 12.566 # length as a string
    constants = [ncells, npart, time_length, num_time_steps, length]
    
    # type
    type_name = 'LD' # 'LD' or '2SI'
    
    # variable
    scan_var = 'num_time_steps'#'ncells'/'npart'/'time_length'/'num_time_steps'/'length'
    num_repeats = 10
    start_var = 5
    stop_var =  250
    N = 20
    is_log = False
    this_log_base = 10.0
    
    # general 
    is_plot_graphs = False
    is_plot_summary_graph = False
    results = analyse_scan_range(scan_var, num_repeats, start_var, stop_var, N,
                       is_log, constants, type_name,
                       is_plot_graphs, is_plot_summary_graph, this_log_base)
    
    # Fit
    fit_vals, fit_square_errs = perform_fit3(results[0], results[3])
    num_fit_points = 100
    fit_xs = numpy.linspace(start_var, stop_var, num_fit_points)
    fit_ys = pow_func3(fit_vals, fit_xs)
    fit_label = r'Power Fit: $(%.2g\pm %.1g)(x)^{(%.2g\pm %.1g)}+(%.2g\pm %.1g)$' % (
            fit_vals[0],
            numpy.sqrt(fit_square_errs[0]),
            fit_vals[1],
            numpy.sqrt(fit_square_errs[1]),
            fit_vals[2],
            numpy.sqrt(fit_square_errs[2]))
    
    plt.plot(fit_xs, fit_ys, label=fit_label)
    this_label = hf.generate_filename(
            str(type_name), str(length), 'VAR', str(npart),
            str(time_length), str(num_time_steps), 'AVG:{}'.format(num_repeats))
    plt.errorbar(results[0], results[3], yerr=results[4], marker='x', label = this_label)
    plt.ylabel('Time Taken [s]')
    plt.xlabel('Number of Time Steps')
    plt.ioff() # This so that the windows stay open
    plt.legend()
    plt.yscale('linear')
    plt.xscale('linear')
    plt.show()

def plot_time_vs_length():
    # constants
    ncells = 20
    npart = 1000
    time_length = 20
    num_time_steps = 50
    length = 12.566 # length as a string
    constants = [ncells, npart, time_length, num_time_steps, length]
    
    # type
    type_name = 'LD' # 'LD' or '2SI'
    
    # variable
    scan_var = 'length'#'ncells'/'npart'/'time_length'/'num_time_steps'/'length'
    num_repeats = 10
    start_var = 0.5*numpy.pi
    stop_var =  10*numpy.pi
    N = 20
    is_log = False
    this_log_base = 10.0
    
    # general 
    is_plot_graphs = False
    is_plot_summary_graph = False
    results = analyse_scan_range(scan_var, num_repeats, start_var, stop_var, N,
                       is_log, constants, type_name,
                       is_plot_graphs, is_plot_summary_graph, this_log_base)
    
    
    # Fit
    fit_vals, fit_square_errs = perform_fit3(results[0], results[3], variables=[0,-1,0])
    num_fit_points = 100
    fit_xs = numpy.linspace(start_var, stop_var, num_fit_points)
    fit_ys = pow_func3(fit_vals, fit_xs)
    fit_label = r'Power Fit: $(%.2g\pm %.1g)(x)^{(%.2g\pm %.1g)}+(%.2g\pm %.1g)$' % (
            fit_vals[0],
            numpy.sqrt(fit_square_errs[0]),
            fit_vals[1],
            numpy.sqrt(fit_square_errs[1]),
            fit_vals[2],
            numpy.sqrt(fit_square_errs[2]))
    
    plt.plot(fit_xs, fit_ys, label=fit_label)
    
    
    this_label = hf.generate_filename(
            str(type_name), str(length), 'VAR', str(npart),
            str(time_length), str(num_time_steps), 'AVG:{}'.format(num_repeats))
    plt.errorbar(results[0], results[3], yerr=results[4], marker='x', label = this_label)
    plt.ylabel('Time Taken [s]')
    plt.xlabel('Length')
    plt.ioff() # This so that the windows stay open
    plt.legend()
    plt.yscale('linear')
    plt.xscale('linear')
    plt.show()


############################### DAMPING RATES

def plot_damping_vs_length():
    # constants
    ncells = 20
    npart = 1000
    time_length = 20
    num_time_steps = 50
    length = 12.566 # length as a string
    constants = [ncells, npart, time_length, num_time_steps, length]
    
    # type
    type_name = 'LD' # 'LD' or '2SI'
    
    # variable
    scan_var = 'length'#'ncells'/'npart'/'time_length'/'num_time_steps'/'length'
    num_repeats = 1
    start_var = 4*numpy.pi
    stop_var =  4*numpy.pi
    N = 1
    is_log = False
    this_log_base = 10.0
    
    # general 
    is_plot_graphs = True
    is_plot_summary_graph = True
    results = analyse_scan_range(scan_var, num_repeats, start_var, stop_var, N,
                       is_log, constants, type_name,
                       is_plot_graphs, is_plot_summary_graph, this_log_base)
    
    
    # Fit
    fit_vals, fit_square_errs = perform_fit3(results[0], results[3], variables=[0,-1,0])
    num_fit_points = 100
    fit_xs = numpy.linspace(start_var, stop_var, num_fit_points)
    fit_ys = pow_func3(fit_vals, fit_xs)
    fit_label = r'Power Fit: $(%.2g\pm %.1g)(x)^{(%.2g\pm %.1g)}+(%.2g\pm %.1g)$' % (
            fit_vals[0],
            numpy.sqrt(fit_square_errs[0]),
            fit_vals[1],
            numpy.sqrt(fit_square_errs[1]),
            fit_vals[2],
            numpy.sqrt(fit_square_errs[2]))
    
    plt.plot(fit_xs, fit_ys, label=fit_label)
    
    
    this_label = hf.generate_filename(
            str(type_name), str(length), 'VAR', str(npart),
            str(time_length), str(num_time_steps), 'AVG:{}'.format(num_repeats))
    plt.errorbar(results[0], results[3], yerr=results[4], marker='x', label = this_label)
    plt.ylabel('Time Taken [s]')
    plt.xlabel('Length')
    plt.ioff() # This so that the windows stay open
    plt.legend()
    plt.yscale('linear')
    plt.xscale('linear')
    plt.show()

if __name__ == "__main__":
    #plot_scan_noise_vs_numCells()
    #plot_scan_noise_vs_numPart()
    #plot_scan_noise_vs_timeLength()
    #plot_scan_noise_vs_numTimeSteps()
    #plot_scan_noise_vs_length()
    #plot_time_vs_numCells()
    #plot_time_vs_numParts()
    #plot_time_vs_timeLength()
    #plot_time_vs_numTimeSteps()
    #plot_time_vs_length()
    plot_damping_vs_length()
