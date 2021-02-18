# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 14:37:49 2021

@author: dmg530
"""

def find_filename(this_data_run):
    type_name = this_data_run.get_type()
    length = '_L-' +  str(round(this_data_run.get_L(), 3))#round to 3 dp
    num_cells = '_C-' + str(this_data_run.get_ncells())
    num_parts = '_P-' + str(this_data_run.get_npart())
    time_length = '_T-' + str(this_data_run.get_time_length())
    time_steps = '_Tn-' + str(this_data_run.get_num_time_steps())
    iteration_num = '_I-' + str(this_data_run.get_iteration_num())
    name = type_name+length+num_cells+num_parts+time_length+time_steps+iteration_num
    return name

def generate_filename(type_name, L, C, P, T, Tn, I):
    type_name = type_name
    length = '_L-' +  L
    num_cells = '_C-' + C
    num_parts = '_P-' + P
    time_length = '_T-' + T
    time_steps = '_Tn-' + Tn
    iteration_num= '_I-' + I
    name = type_name+length+num_cells+num_parts+time_length+time_steps+iteration_num
    return name