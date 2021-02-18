#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 13:54:14 2021

@author: dmg530
"""

import numpy 
data = 10*numpy.random.random(100)
data.sort()
bins = numpy.linspace(0, 10, 10)
digitized = numpy.digitize(data, bins)
print(data)
print(digitized)
offset = data-digitized+1
print(offset)

a = numpy.zeros([1])
print(a.type)