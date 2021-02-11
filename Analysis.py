# -*- coding: utf-8 -*-
"""
Class for an Analysis of a Data Run
"""
from Summary import Summary

from numpy import arange, linspace
from numpy import argmin, argmax

from scipy import interpolate
import matplotlib.pyplot as plt # Matplotlib plotting library
class Analysis:
    """
    Displays three plots: phase space, charge density, and velocity distribution
    """
    def __init__(self, this_Summary):
        s = this_Summary
        # Summary stores an array of the first-harmonic amplitude
        self.analyse_first_harmonic_time(s.t, s.firstharmonic)
        self.analyse_first_harmonic_time(s.t, s.firstharmonic_no_abs)
        # Make a semilog plot to see exponential damping
        self.plot_semilog(s.t, s.firstharmonic)
    
    def plot_semilog(self, x, y):
        plt.figure()
        plt.plot(x,y)
        plt.xlabel("Time [Normalised]")
        plt.ylabel("First harmonic amplitude [Normalised]")
        plt.yscale('log')
        
        plt.ioff() # This so that the windows stay open
        plt.show()
        
    
    def is_pos_diff(self, vals, test_index):
        is_this_pos_grad = None
        diff = vals[test_index+1]-vals[test_index]
        if diff >= 0:
            is_this_pos_grad = True
        else:
            is_this_pos_grad = False
        return is_this_pos_grad
    
    def turn_point(self, xs, interp, target_index, is_min):
        N = 10
        LH = xs[target_index-1]
        RH = xs[target_index+ 1]
        x = linspace(LH, RH, num=N)
        y = interp(x)
        turn_index = None
        if is_min:
            turn_index = argmin(y)
        else:
            turn_index = argmax(y)
        turn_point_y = y[turn_index]
        turn_point_x = x[turn_index]
        return turn_point_x, turn_point_y
    
    # for finding minimum and maximum of first harmonic - time spectrum
    def analyse_first_harmonic_time(self, ts, fhs):
        max_pos_indexes = [1] # the first is a peak
        min_pos_indexes = []
        is_prev_pos_grad = self.is_pos_diff(fhs, 0)
        for i in range(len(fhs)-1):
            # check if this one is a positive gradient
            is_this_pos_grad = self.is_pos_diff(fhs, i)
            # check for change from previous
            if is_prev_pos_grad and not is_this_pos_grad:
                max_pos_indexes.append(i)
            if not is_prev_pos_grad and is_this_pos_grad:
                min_pos_indexes.append(i)
            is_prev_pos_grad = is_this_pos_grad
        f = interpolate.interp1d(ts, fhs, kind='quadratic')
        x_minima, y_minima, x_maxima, y_maxima = [], [], [], []
        for min_index in min_pos_indexes:
            x,y = self.turn_point(ts, f, min_index, True)
            x_minima.append(x)
            y_minima.append(y)
            plt.plot(x, y, 'x', color='r')
        for max_index in max_pos_indexes:
            x,y = self.turn_point(ts, f, max_index, False)
            x_maxima.append(x)
            y_maxima.append(y)
            plt.plot(x, y, 'o', color='k')
        
        xnew = arange(ts[0], ts[-1], 0.05)
        ynew = f(xnew)   # use interpolation function returned by `interp1d`
        plt.plot(ts, fhs, xnew, ynew, '-', label='Curve Fit in Peak Region')
        plt.plot(ts, fhs, label='Underlying Data')
        plt.xlabel("Time [Normalised]")
        plt.ylabel("First harmonic amplitude [Normalised]")
        plt.legend()
        plt.ioff() # This so that the windows stay open
        plt.show()
        