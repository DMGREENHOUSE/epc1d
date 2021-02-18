#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 00:43:46 2021

@author: dmg530
"""

#from numpy import zeros, ndarray
import numpy as np
cimport numpy as np

def calc_density_old(position, ncells, L):
    """ Calculate charge density given particle positions
    
    Input
      position  - Array of positions, one for each particle
                  assumed to be between 0 and L
      ncells    - Number of cells
      L         - Length of the domain

    Output
      density   - contains 1 if evenly distributed
    """
    cdef double dx = L / ncells       # Uniform cell spacing
    
    cdef int nparticles = len(position)
    cdef np.ndarray[double, ndim=1] density = np.zeros([ncells])
    
    for p in position / dx:    # Loop over all the particles, converting position into a cell number
        plower = int(p)        # Cell to the left (rounding down)
        offset = p - plower    # Offset from the left
        density[plower] += 1. - offset
        density[(plower + 1) % ncells] += offset
    # nparticles now distributed amongst ncells
    density *= float(ncells) / float(nparticles)  # Make average density equal to 1
    return density





















    """
    for (i = 0; i < nparticles_c; ++i)
    {
        p = position_cell_c[i]
        plower = (int)floor(p)        # Cell to the left (rounding down)
        offset = p - plower    # Offset from the left
        density_c[plower] += 1. - offset
        density_next_c[plower+1] += offset
    }
    #density_next = numpy.concatenate(([density_next[-1]],density_next[0:-1]))
    #density += density_next[0:-1]
    #density[0] += density_next[-1]
    # nparticles now distributed amongst ncells
    cdef float norm = ncells_c / nparticles_c
    density_c *= norm  # Make average density equal to 1
    
    return density_c
 """"""
    cdef numpy.ndarray density = numpy.zeros([ncells])
    cdef numpy.ndarray density_next = numpy.zeros([ncells+1])
    cdef int nparticles = len(position)
    
    cdef int ncells_c = ncells
    cdef float [:] position_c = position
    cdef float density_c[20]
    
    cdef float dx = L / ncells       # Uniform cell spacing
    
    cdef float p, offset
    cdef int plower
    
    for p in position / dx:    # Loop over all the particles, converting position into a cell number
        plower = int(p)        # Cell to the left (rounding down)
        offset = p - plower    # Offset from the left
        density[plower] += 1. - offset
        density_next[plower+1] += offset
    density_next = numpy.concatenate(([density_next[-1]],density_next[0:-1]))
    density += density_next[0:-1]
    density[0] += density_next[-1]
    # nparticles now distributed amongst ncells
    density *= ncells / nparticles  # Make average density equal to 1
    
    return density
    """