# -*- coding: utf-8 -*-
"""
Class for an Analysis of a Data Run
"""
from Summary import Summary

from numpy import arange, linspace
from numpy import argmin, argmax, exp, zeros

from scipy import interpolate, optimize
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
        
    def exponential_func(self, fit_vals,t):
        return fit_vals[0]*exp(-fit_vals[1]*t)+fit_vals[2]
    
    def check_covariances(self,cov_vals):
        is_above_one = False
        check_val = 1
        for i in cov_vals:
            for j in i:
                if j>check_val:
                    is_above_one = True
        return is_above_one
    
    def perform_fit(self,func, xmaxs,ymaxs,threes):
        fit_vals, cov_vals = optimize.curve_fit(
                lambda t,a,b,c: func([a,b,c],t),
                xmaxs,  ymaxs,  p0=(threes[0], threes[1],threes[2]))
        return fit_vals, cov_vals
    
    def fit(self,is_max, xmaxs,ymaxs):
        func = self.exponential_func
        xmaxs,ymaxs = self.find_suitable_points(is_max, xmaxs, ymaxs)
        num_found = len(xmaxs)
        SKIP = False
        if num_found <= 3:
            print('Not enough points found')
            SKIP = True
        else:
            pass
        count = 0
        count_lim = 3
        nulls = [[1,1,1],
                 [-1,1,0],
                 [-1,-1,-1]
                ]
        fit_vals = None
        if not SKIP:
            while count<count_lim:
                # try with first fit
                fit_vals, cov_vals = self.perform_fit(func, xmaxs,ymaxs,nulls[count])
                print(cov_vals)
                is_bad_fit = False
                if not SKIP:
                    is_bad_fit = self.check_covariances(cov_vals)
                print(is_bad_fit)
                if is_bad_fit:
                    count+=1
                else:
                    count=4
            print('----COMPLETE---')
        if count == 3 or SKIP:
            # no good fit was found
            print('None found', fit_vals)
            return None
        else:
            # a good fit was found
            return fit_vals
    
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
    
    def fill_array(self, func, fit, x):
        N = len(x)
        y = zeros([N])
        for i in range(N):
            y[i] = func(fit, x[i])
        
        return y
    
    def find_suitable_points(self, is_max, initial_x_points, initial_y_points):
        # we assume the first two points are true
        diff = initial_x_points[1] - initial_x_points[0] 
        BUFFER_PERCENTAGE = 0.25 # allow some buffer
        new_x_points = initial_x_points[0:1]
        print('!!!, ', new_x_points)
        new_y_points = initial_y_points[0:1]
        num_points = len(initial_x_points)
        for i in range(num_points-1):
            this_point = initial_x_points[i+1]
            this_diff = this_point-new_x_points[-1] #compare to the final 'good' point
            print(this_point, new_x_points[-1], this_diff)
            if this_diff < (1+BUFFER_PERCENTAGE)*diff and this_diff > (1-BUFFER_PERCENTAGE)*diff:
                potential_y = initial_y_points[i+1]
                if is_max and potential_y < new_y_points[-1]:
                    new_x_points.append(this_point)
                    new_y_points.append(potential_y)
                elif not is_max and potential_y > new_y_points[-1]:
                    new_x_points.append(this_point)
                    new_y_points.append(potential_y)
            else:
                print('woo')
        print(initial_x_points)
        print(new_x_points)
        plt.plot(new_x_points, new_y_points, 'o', color='k', label='Suitables')
        return new_x_points, new_y_points
        
    
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
            #plt.plot(x, y, 'x', color='r')
        for max_index in max_pos_indexes:
            x,y = self.turn_point(ts, f, max_index, False)
            x_maxima.append(x)
            y_maxima.append(y)
        
        plt.plot(x_minima, y_minima, 'x', color='r', label='Minima')
        plt.plot(x_maxima, y_maxima, '.', color='b', label='Maxima')
        minima_fits = self.fit(False, x_minima, y_minima)
        maxima_fits = self.fit(True, x_maxima, y_maxima)
        
        xnew = arange(ts[0], ts[-1], 0.05)
        ynew = f(xnew)   # use interpolation function returned by `interp1d`
        
        if minima_fits is not None:
            y_min_fits = self.fill_array(self.exponential_func,minima_fits, ts)
            plt.plot(ts, y_min_fits)
        if maxima_fits is not None:
            y_max_fits = self.fill_array(self.exponential_func,maxima_fits, ts)
            plt.plot(ts, y_max_fits)
        plt.plot(ts, fhs, xnew, ynew, '-', label='Curve Fit in Peak Region')
        plt.plot(ts, fhs, label='Underlying Data')
        plt.xlabel("Time [Normalised]")
        plt.ylabel("First harmonic amplitude [Normalised]")
        plt.legend()
        plt.ioff() # This so that the windows stay open
        plt.show()
        