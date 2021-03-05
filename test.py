#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 13:42:30 2021

@author: dmg530
"""
import numpy as np
density_c = np.array([0,1,2,3,4,5,6,7,8,9])
density_next_c = np.array([0,0,10,20,30,40,50,60,70,80,90])
density_c += density_next_c[0:-1]
print(density_c)
density_c[0] += density_next_c[-1]

print(density_c)