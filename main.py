# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 13:49:06 2021

@author: dmg530
"""

from Data_Run import Data_Run
from Analysis import Analysis
import pickle
from helper_functions import *
import time

start = time.time()
new_data_run = Data_Run(ncells=20, npart=2000)
end = time.time()
time_taken = end-start
print('TIME TAKEN: ', time_taken)
this_summary = new_data_run.get_summary()
#new_analysis = Analysis(this_summary)
"""
filename = find_filename(new_data_run)
filepath = r'Data/'+filename


outfile = open(filepath, 'wb')
pickle.dump(this_summary,outfile)
outfile.close()

# input the object
infile = open(filename,'rb')
new_summary = pickle.load(infile)
infile.close()
new_analysis = Analysis(new_summary)
"""