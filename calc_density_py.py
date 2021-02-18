# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 13:53:35 2021

@author: dmg530
"""
from numpy import zeros, array, concatenate, arange, floor, ones
from pandas import DataFrame, cut

def calc_density(position, ncells, L):
    """ Calculate charge density given particle positions
    
    Input
      position  - Array of positions, one for each particle
                  assumed to be between 0 and L
      ncells    - Number of cells
      L         - Length of the domain

    Output
      density   - contains 1 if evenly distributed
    """
    # This is a crude method and could be made more efficient
    density = zeros([ncells])
    nparticles = len(position)
    
    dx = L / ncells       # Uniform cell spacing
    OLD = True
    #if (ncells>nparticles):
    if OLD:
    
        for p in position / dx:    # Loop over all the particles, converting position into a cell number
            plower = int(p)        # Cell to the left (rounding down)
            offset = p - plower    # Offset from the left
            density[plower] += 1. - offset
            density[(plower + 1) % ncells] += offset
        # nparticles now distributed amongst ncells
        density *= float(ncells) / float(nparticles)  # Make average density equal to 1
        """
    else:
        density_next = zeros([ncells+1])
        for p in position / dx:    # Loop over all the particles, converting position into a cell number
            plower = int(p)        # Cell to the left (rounding down)
            offset = p - plower    # Offset from the left
            density[plower] += 1. - offset
            density_next[plower+1] += offset
        density_next = concatenate(([density_next[-1]],density_next[0:-1]))
        density += density_next[0:-1]
        density[0] += density_next[-1]
        # nparticles now distributed amongst ncells
        density *= float(ncells) / float(nparticles)  # Make average density equal to 1
        """
    else:
        cells = arange(ncells+1)
        position/=dx # convert position into a cell number
        counts = ones(nparticles)
        df = DataFrame({"positions": position, "count":counts})
        df["bin"] = cut(df["positions"], cells) # bin into the cells
        df = df.groupby('bin').sum() # sum up bins
        df = df.sort_values(by="bin")
        summed = df["positions"].values
        count = df["count"].values
        base_column = (cells[0:-1]+1)*count-summed
        next_column = concatenate(([summed[-1]],summed[0:-1]))
        
        density = base_column+next_column
        density *= float(ncells) / float(nparticles)  # Make average density equal to 1
    
    """
    else:
        cells = arange(ncells+1)
        position/=dx # convert position into a cell number
        df = DataFrame({"positions": position})
        offset = df['positions'] - df['positions'].apply(floor)
        df['this']= -offset+1 # the amount in this bin
        df['next']=offset # the amount in the next bin
        df["bin"] = cut(df["positions"], cells) # bin into the cells
        df = df.sort_values(by="bin")
        df = df.groupby('bin').sum() # sum up bins
        base_column = df["this"].values # particles in these cells
        next_column = df["next"].values # particles in next cells
        # shift next col to align with this col and move last val to first
        next_column = concatenate(([next_column[-1]],next_column[0:-1]))
        density = base_column+next_column
        density *= float(ncells) / float(nparticles)  # Make average density equal to 1
    """
    """
        density = zeros([ncells+1])
        position/=dx # convert position into a cell number
        for p in position:    # Loop over all the particles
            plower = int(p)        # Cell to the left (rounding down)
            offset = p - plower    # Offset from the left
            density[plower] += 1. - offset
            density[(plower + 1)] += offset
        density[-2]+=density[-1]
        density = density[0:-1]
        # nparticles now distributed amongst ncells
        density *= float(ncells) / float(nparticles)  # Make average density equal to 1
        """
    """
        # arange the particles in order
        position.sort()
        position/=dx #convert distance to cell number
        for cell in range(ncells):
            count=0
            max_left = len(position) # the maximum number of remaining particles
            while (count<max_left) and position[count]<(cell+1): # only look at the start of this ordered list
                offset = position[count] - (cell)
                density[cell] += 1. - offset
                density[(cell + 1) % ncells] += offset
                count+=1
            position = position[count:]# shrink position array since they've already been counted
        density *= float(ncells) / float(nparticles)  # Make average density equal to 1
        """
    return density