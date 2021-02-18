# -*- coding: utf-8 -*-
"""
Class for an Analysis of a Data Run
"""
from Summary import Summary

from numpy import arange, linspace, log, mean
from numpy import argmin, argmax, exp, zeros

from scipy import interpolate, optimize
import matplotlib.pyplot as plt # Matplotlib plotting library
class Analysis:
    """
    Displays three plots: phase space, charge density, and velocity distribution
    """
    def __init__(self, this_Summary, is_plot_graphs):
        s = this_Summary
        self.time_taken = s.get_time_taken()
        is_print_results = False
        is_print_data_point_warning = False
        
        # Summary stores an array of the first-harmonic amplitude
        
        # Find Analysis for semi-log approach
        self.log_noise_floor, self.log_damping_rate = self.analyse_first_harmonic_time(
                s.t, s.firstharmonic, True, True, is_plot_graphs, is_print_data_point_warning)
        
        # Find Analysis for no log approach
        self.nolog_noise_floor, self.nolog_damping_rate  = self.analyse_first_harmonic_time(
                s.t, s.firstharmonic, True, False, is_plot_graphs, is_print_data_point_warning)
        
        # Find Analysis for no log, no abs approach
        self.noabs_noise_floor, self.noabs_damping_rate  = self.analyse_first_harmonic_time(
                s.t, s.firstharmonic_no_abs, False, False, is_plot_graphs, is_print_data_point_warning)
        
        if is_print_results:
            print('Noise level according to log: ', self.log_noise_floor)
            print('Damping Rate according to log: ', self.log_damping_rate)
            print('Noise level according to combined no log: ', self.nolog_noise_floor)
            print('Damping Rate according to combined no log: ', self.nolog_damping_rate)
            print('Noise level according to max/min separated no log: ', self.noabs_noise_floor)
            print('Damping Rate according to max/min separated no log: ', self.noabs_damping_rate)
        
        # Make a semilog plot to see exponential damping
        if is_plot_graphs:
            self.plot_semilog(s.t, s.firstharmonic)
    
    def get_results_time(self):
        return self.time_taken
    
    def get_results_log(self):
        return [self.log_noise_floor, self.log_damping_rate]
    
    def get_results_nolog(self):
        return [self.nolog_noise_floor, self.nolog_damping_rate]
    
    def get_results_noabs(self):
        return [self.noabs_noise_floor, self.noabs_damping_rate]
    
    def exponential_linear_func(self, fit_vals,t):
        #note linear under log
        return fit_vals[0]*exp(-fit_vals[1]*t)
        
    def exponential_func(self, fit_vals,t):
        return fit_vals[0]*exp(-fit_vals[1]*t)+fit_vals[2]
    
    def check_covariances(self,cov_vals):
        if cov_vals is None:
            return True
        else:
            is_above_one = False
            check_val = 1
            for i in cov_vals:
                for j in i:
                    if j>check_val:
                        is_above_one = True
            return is_above_one
    
    def perform_two_fit(self,func, xmaxs,ymaxs,twos):
        try:
            fit_vals, cov_vals = optimize.curve_fit(
                    lambda t,a,b: func([a,b],t),
                    xmaxs,  ymaxs,  p0=(twos[0], twos[1]), maxfev=1000)
            return fit_vals, cov_vals
        except RuntimeError:
            return None, None
    
    def perform_three_fit(self,func, xmaxs,ymaxs,threes):
        try:
            fit_vals, cov_vals = optimize.curve_fit(
                    lambda t,a,b,c: func([a,b,c],t),
                    xmaxs,  ymaxs,  p0=(threes[0], threes[1],threes[2]), maxfev=1000)
            return fit_vals, cov_vals
        except RuntimeError:
            return None, None
    
    def fit(self, is_semi_log, func, is_max, xmaxs,ymaxs, is_plot_graphs, is_print_data_point_warning):
        nulls_two = [[-2,1],
                 [0,1],
                 [1,-1]
                ]
        nulls_three = [[2,1,1],
                 [2,-1,0],
                 [-2,-1,-1]
                ]
            
        xmaxs,ymaxs, average_noise = self.find_suitable_points(is_max, xmaxs, ymaxs, is_plot_graphs)
        num_found = len(xmaxs)
        SKIP = False
        if num_found <= 3:
            if is_print_data_point_warning:
                print('Not enough points found')
            SKIP = True
        else:
            pass
        count = 0
        count_lim = 3
        fit_vals = None
        cov_vals = None
        if not SKIP:
            while count<count_lim:
                # try with first fit
                if is_semi_log:
                    fit_vals, cov_vals = self.perform_two_fit(func, xmaxs,ymaxs,nulls_two[count])
                else:
                    fit_vals, cov_vals = self.perform_three_fit(func, xmaxs,ymaxs,nulls_three[count])
                is_bad_fit = False
                if not SKIP:
                    is_bad_fit = self.check_covariances(cov_vals)
                #print(is_bad_fit)
                if is_bad_fit:
                    count+=1
                else:
                    count=4
        if count == 3 or SKIP:
            # no good fit was found
            #print('None found', fit_vals)
            return None, average_noise
        else:
            # a good fit was found
            return fit_vals, average_noise
    
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
        
    def find_suitable_points(self, is_max, initial_x_points, initial_y_points, is_plot_graphs):
        # we assume the first two points are true
        diff = 0
        if len(initial_x_points)>1:
            diff = initial_x_points[1] - initial_x_points[0] 
        BUFFER_PERCENTAGE = 0.15 # allow some buffer
        suitable_x_points = initial_x_points[0:1]
        suitable_y_points = initial_y_points[0:1]
        noise_x_points = []
        noise_y_points = []
        
        num_points = len(initial_x_points)
        for i in range(num_points-1):
            this_point = initial_x_points[i+1]
            potential_y = initial_y_points[i+1]
            this_diff = this_point-suitable_x_points[-1] #compare to the final 'good' point
            if this_diff < (1+BUFFER_PERCENTAGE)*diff and this_diff > (1-BUFFER_PERCENTAGE)*diff:
                if is_max and potential_y < suitable_y_points[-1]:
                    suitable_x_points.append(this_point)
                    suitable_y_points.append(potential_y)
                elif not is_max and potential_y > suitable_y_points[-1]:
                    suitable_x_points.append(this_point)
                    suitable_y_points.append(potential_y)
            else:
                noise_x_points.append(this_point)
                noise_y_points.append(potential_y)
        average_noise = mean(noise_y_points)
        if is_plot_graphs and is_max:         
            plt.plot(suitable_x_points, suitable_y_points, 'o', color='k', label='Suitable Maxima')
            plt.plot(noise_x_points, noise_y_points, 'x', color='k', label='Noise Maxima')
            plt.axhline(y=average_noise, linestyle='--', color='k', label='Average Noise Maxima')
        elif is_plot_graphs and not is_max:
            plt.plot(suitable_x_points, suitable_y_points, 'o', color='r', label='Suitable Minima')
            plt.plot(noise_x_points, noise_y_points, 'x', color='r', label='Noise Minima')
            plt.axhline(y=average_noise, linestyle='--', color='r', label='Average Noise Maxima')
            
        return suitable_x_points, suitable_y_points, average_noise
        
    
    # for finding minimum and maximum of first harmonic - time spectrum
    def analyse_first_harmonic_time(self, ts, fhs, is_abs, is_semi_log, is_plot_graphs, is_print_data_point_warning):
        func = None
        if is_semi_log:
            plt.yscale('log')
            func = self.exponential_linear_func
        else:
            func = self.exponential_func
            
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
        if not is_abs:
            for min_index in min_pos_indexes:
                x,y = self.turn_point(ts, f, min_index, True)
                x_minima.append(x)
                y_minima.append(y)
                #plt.plot(x, y, 'x', color='r')
        for max_index in max_pos_indexes:
            x,y = self.turn_point(ts, f, max_index, False)
            x_maxima.append(x)
            y_maxima.append(y)
        
        #plt.plot(x_minima, y_minima, 'x', color='r', label='Minima')
        #plt.plot(x_maxima, y_maxima, '.', color='b', label='Maxima')
        
        
        minima_fits = []
        minima_average_noise = None
        if not is_abs:
            minima_fits, minima_average_noise = self.fit(is_semi_log, func, False, x_minima, y_minima, is_plot_graphs, is_print_data_point_warning)
        maxima_fits, maxima_average_noise = self.fit(is_semi_log, func, True, x_maxima, y_maxima, is_plot_graphs, is_print_data_point_warning)
        
        noise_floor = self.find_noise_floor(is_abs, is_semi_log,
                              minima_fits, minima_average_noise,
                              maxima_fits, maxima_average_noise)
        damping_rate = self.find_damping_rate(is_abs, is_semi_log,
                              minima_fits, minima_average_noise,
                              maxima_fits, maxima_average_noise)
        
        xnew = arange(ts[0], ts[-1], 0.05)
        ynew = f(xnew)   # use interpolation function returned by `interp1d`
        
        
        if is_plot_graphs:
            if not is_abs:
                if minima_fits is not None:
                    y_min_fits = self.fill_array(func,minima_fits, ts)
                    plt.plot(ts, y_min_fits, label='Fit of Suitable Minima Points')
            if maxima_fits is not None:
                y_max_fits = self.fill_array(func,maxima_fits, ts)
                plt.plot(ts, y_max_fits, label='Fit of Suitable Maxima Points')
            plt.plot(xnew, ynew, '-', label='Curve Fit in Peak Region')
            plt.plot(ts, fhs, label='Underlying Data')
            plt.xlabel("Time [Normalised]")
            plt.ylabel("First harmonic amplitude [Normalised]")
            plt.legend()
            plt.ioff() # This so that the windows stay open
            plt.show()
        return noise_floor, damping_rate
        
    def find_noise_floor(self, is_abs, is_semi_log,
                              minima_fits, minima_average_noise,
                              maxima_fits, maxima_average_noise):
        noise_floor = None
        if is_semi_log:
            # here we want to just take the noise maxima average
            noise_floor = maxima_average_noise
        else:
            if is_abs:
                # here we want to pick between both the exponential floor and maxima average
                if maxima_fits is None:
                    noise_floor = maxima_average_noise
                else:
                    noise_floor = max(maxima_fits[2], maxima_average_noise)
            elif (not is_abs):
                # here we want to pick between both the exponential floor and minima average
                #     for both the maxima and minima and then take the average
                noise_for_maxima = None
                noise_for_minima = None
                if (maxima_fits is None):
                    noise_for_maxima = maxima_average_noise
                else:
                    noise_for_maxima = max(maxima_fits[2], maxima_average_noise)
                if (minima_fits is None):
                    noise_for_minima = minima_average_noise
                else:
                    noise_for_minima = min(minima_fits[2], minima_average_noise)
                    
                diff_max_min = abs(noise_for_maxima - abs(noise_for_minima))
                noise_floor = (noise_for_maxima+abs(noise_for_minima))/2
                BUFFER = 1
                if diff_max_min/noise_floor > BUFFER:
                    print('Difference between max and min noise floors for non abs exceeds buffer.\n \t Difference is: ',
                          diff_max_min)
                    
        return noise_floor
        
    def find_damping_rate(self, is_abs, is_semi_log,
                              minima_fits, minima_average_noise,
                              maxima_fits, maxima_average_noise):
        damping_rate = None
        if is_semi_log and (maxima_fits is not None):
            # here the rate is simply the gradient of the slope
            damping_rate = maxima_fits[1]
        else:
            if is_abs and (maxima_fits is not None):
                # here we want B in y=Ae^(Bx)+C
                damping_rate = maxima_fits[1]
            elif (not is_abs) and (maxima_fits is not None) and  (minima_fits is not None):
                # here we want B in y=Ae^(Bx)+C
                #     for both the maxima and minima and then take the average
                damping_rate_for_maxima = maxima_fits[1]
                damping_rate_for_minima = minima_fits[1]
                diff_max_min = abs(abs(damping_rate_for_maxima) - abs(damping_rate_for_minima))
                damping_rate = (damping_rate_for_maxima-damping_rate_for_minima)/2
                BUFFER = 0.1
                if diff_max_min/damping_rate > BUFFER:
                    print('Difference between max and min damping rates for non abs exceeds buffer.\n \t Difference is: ',
                          diff_max_min)
        return damping_rate