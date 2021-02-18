#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 15:03:30 2021

@author: dmg530
"""

from setuptools import setup
from Cython.Build import cythonize

setup(
      name = 'Cython Modules',
      ext_modules = cythonize(["*.pyx"]),# = cythonize("calc_density.pyx", "calc_density_adapted.pyx", "calc_density_old.pyx")
)