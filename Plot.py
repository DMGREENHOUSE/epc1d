# -*- coding: utf-8 -*-
"""
Plot Class
"""

from numpy import linspace
from numpy import sqrt, histogram

import matplotlib.pyplot as plt # Matplotlib plotting library

try:
    import matplotlib.gridspec as gridspec  # For plot layout grid
    got_gridspec = True
except:
    got_gridspec = False

from calc_density import calc_density

class Plot:
    """
    Displays three plots: phase space, charge density, and velocity distribution
    """
    def __init__(self, pos, vel, ncells, L,is_plot_animation):
        
        d = calc_density(pos, ncells, L)
        vhist, bins  = histogram(vel, int(sqrt(len(vel))))
        vbins = 0.5*(bins[1:]+bins[:-1])
        self.is_plot_animation = is_plot_animation
        # Plot initial positions
        if self.is_plot_animation:
            if got_gridspec:
                self.fig = plt.figure()
                self.gs = gridspec.GridSpec(4, 4)
                ax = self.fig.add_subplot(self.gs[0:3,0:3])
                self.phase_plot = ax.plot(pos, vel, '.')[0]
                ax.set_title("Phase space")
                
                ax = self.fig.add_subplot(self.gs[3,0:3])
                self.density_plot = ax.plot(linspace(0, L, ncells), d)[0]
                
                ax = self.fig.add_subplot(self.gs[0:3,3])
                self.vel_plot = ax.plot(vhist, vbins)[0]
            else:
                self.fig = plt.figure()
                self.phase_plot = plt.plot(pos, vel, '.')[0]
                
                self.fig = plt.figure()
                self.density_plot = plt.plot(linspace(0, L, ncells), d)[0]
                
                self.fig = plt.figure()
                self.vel_plot = plt.plot(vhist, vbins)[0]
            plt.ion()
            plt.show()
        
    def __call__(self, pos, vel, ncells, L, t):
        if self.is_plot_animation:
            d = calc_density(pos, ncells, L)
            vhist, bins  = histogram(vel, int(sqrt(len(vel))))
            vbins = 0.5*(bins[1:]+bins[:-1])
            
            self.phase_plot.set_data(pos, vel) # Update the plot
            self.density_plot.set_data(linspace(0, L, ncells), d)
            self.vel_plot.set_data(vhist, vbins)
            plt.draw()
            plt.pause(0.05)
        else:
            pass