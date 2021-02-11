# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 13:49:06 2021

@author: dmg530
"""

from Data_Run import Data_Run
from Analysis import Analysis

new_data_run = Data_Run()
this_summary = new_data_run.get_summary()
new_analysis = Analysis(this_summary)