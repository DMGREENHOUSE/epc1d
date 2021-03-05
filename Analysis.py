# -*- coding: utf-8 -*-
"""
Class for an Analysis of a Data Run
"""
from Summary import Summary

from numpy import arange, linspace, log, mean, pi, concatenate
from numpy import argmin, argmax, exp, zeros, array, std, sort
from math import isnan

from scipy import interpolate, optimize
import matplotlib.pyplot as plt # Matplotlib plotting library
class Analysis:
    """
    Displays three plots: phase space, charge density, and velocity distribution
    """
    def __init__(self, this_Summary, is_plot_graphs, METHOD='Method 2'):
        
        s = this_Summary
        # find the method to be used
        is_abs = None
        is_log = None
        firstharmonic = None
        if METHOD == 'Method 1':
            is_abs = True
            is_log = True
            firstharmonic = s.firstharmonic
        elif METHOD == 'Method 2':
            is_abs = True
            is_log = False
            firstharmonic = s.firstharmonic
        elif METHOD == 'Method 3':
            is_abs = False
            is_log = False
            firstharmonic = s.firstharmonic_no_abs
        else:
            print('ERROR: invalid method')
        # this was a late implementation so only get if possible
        self.time_taken = None
        try:
            self.time_taken = s.get_time_taken()
        except:
            pass
        is_print_results = False
        is_print_data_point_warning = False
        is_plot_semi_log = False
        # Summary stores an array of the first-harmonic amplitude
        
        # Find Analysis for method
        self.noise_floor, self.damping_rate, self.frequency, self.frequency_err = self.analyse_first_harmonic_time(
                s.t, firstharmonic, is_abs, is_log, is_plot_graphs, is_print_data_point_warning)
        
        
        if is_print_results:
            print('Noise level according to {}: '.format(METHOD), self.noise_floor)
            print('Damping Rate according to {}: '.format(METHOD), self.damping_rate)
            print('Angular Frequency according to {}: '.format(METHOD), self.frequency, '+/-', self.frequency_err)
        
        # Make a semilog plot to see exponential damping
        if is_plot_semi_log:
            self.plot_semilog(s.t, s.firstharmonic)

    def get_noise_floor(self):
        return self.noise_floor

    def get_damping_rate(self):
        return self.damping_rate

    def get_frequency(self):
        return self.frequency

    def get_frequency_err(self):
        return self.frequency_err
    
    def get_time(self):
        return self.time_taken
    
    def get_results(self):
        return [self.noise_floor,
                self.damping_rate,
                self.time_taken,
                self.frequency,
                self.frequency_err
                ]
    
    
    
    
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
                 [-2,1,-1]
                ]
            
        suit_xs, suit_ys, average_noise = self.find_suitable_points(is_max, xmaxs, ymaxs, is_plot_graphs)
        num_found = len(suit_xs)
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
                    fit_vals, cov_vals = self.perform_two_fit(func, suit_xs,suit_ys,nulls_two[count])
                else:
                    fit_vals, cov_vals = self.perform_three_fit(func, suit_xs,suit_ys,nulls_three[count])
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
            return None, average_noise, suit_xs
        else:
            # a good fit was found
            return fit_vals, average_noise, suit_xs
    
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
    
    def turn_point_x_range(self, xs, target_index):
        # find range in which to look for a turn point
        N = 50
        LH = xs[target_index-1]
        RH = xs[target_index+ 1]
        if target_index==1:
            #shift search for first point to allow for extrapolation
            LH= - 3*(xs[1]-xs[0])
        xs = linspace(LH, RH, num=N)
        return xs
    
    def turn_point(self, xs, interp, target_index, is_min):
        x = self.turn_point_x_range(xs, target_index)
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
            else:
                noise_x_points.append(this_point)
                noise_y_points.append(potential_y)
        average_noise = mean(noise_y_points)
        if isnan(average_noise):
            average_noise = 0
            
        # was a suitable point but now below the average noise
        
        # remove 'suitables' if below the noise floor
        suit_xs, suit_ys, no_longer_xs, no_longer_ys = self.check_suitables_for_noise_floor(suitable_x_points, suitable_y_points, average_noise, is_max)
        
        
        if is_plot_graphs and is_max:         
            plt.plot(suit_xs, suit_ys, 'o', color='k', label='Suitable Maxima')
            plt.plot(no_longer_xs, no_longer_ys, '^', color='k', label='No Longer Suitable Maxima')
            plt.plot(noise_x_points, noise_y_points, 'x', color='k', label='Noise Maxima')
            plt.axhline(y=average_noise, linestyle='--', color='k', label='Average Noise Maxima')
        elif is_plot_graphs and not is_max:
            plt.plot(suit_xs, suit_ys, 'o', color='r', label='Suitable Minima')
            plt.plot(no_longer_xs, no_longer_ys, '^', color='r', label='No Longer Suitable Minima')
            plt.plot(noise_x_points, noise_y_points, 'x', color='r', label='Noise Minima')
            plt.axhline(y=average_noise, linestyle='--', color='r', label='Average Noise Maxima')
            
        return suit_xs, suit_ys, average_noise
        
    def check_suitables_for_noise_floor(self,suitable_x_points, suitable_y_points, average_noise, is_max):
        new_suitable_xs = []
        new_suitable_ys = []
        no_longer_suitable_xs = []
        no_longer_suitable_ys = []
        if is_max:
            for i in range(len(suitable_y_points)):
                if suitable_y_points[i] >= average_noise:
                    new_suitable_xs.append(suitable_x_points[i])
                    new_suitable_ys.append(suitable_y_points[i])
                else:
                    no_longer_suitable_xs.append(suitable_x_points[i])
                    no_longer_suitable_ys.append(suitable_y_points[i])
        else:
            for i in range(len(suitable_y_points)):
                if suitable_y_points[i] <= average_noise:
                    new_suitable_xs.append(suitable_x_points[i])
                    new_suitable_ys.append(suitable_y_points[i])
                else:
                    no_longer_suitable_xs.append(suitable_x_points[i])
                    no_longer_suitable_ys.append(suitable_y_points[i])
            
        return new_suitable_xs, new_suitable_ys, no_longer_suitable_xs, no_longer_suitable_ys
    
    def coarsley_find_min_max_indices(self, fhs):
        max_pos_indexes = [1] # the first is a peak - use 1 because we will use 0 for fit
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
        return max_pos_indexes, min_pos_indexes
    
    def curve_fit_min_max_regions(self, pos_indexes, ts, fhs):
        local_fits = []
        for i in range(len(pos_indexes)):
            index = pos_indexes[i]
            RANGE = 2 # give a couple of points either side to draw quadratic fit
            # deal with boundaries
            if index-RANGE<0:
                xs = ts[:index+RANGE]
                ys = fhs[:index+RANGE]
            elif index+RANGE>=len(ts):
                xs = ts[index-RANGE:]
                ys = fhs[index-RANGE:]
            else:
                xs = ts[index-RANGE:index+RANGE]
                ys = fhs[index-RANGE:index+RANGE]
            
            try:
                # attempt to allow extrapolation
                local_fit = interpolate.interp1d(xs, ys, kind='quadratic', bounds_error=False,fill_value='extrapolate')
                local_fits.append(local_fit)
            except:
                local_fit = interpolate.interp1d(xs, ys, kind='quadratic')
                local_fits.append(local_fit)
        return local_fits
    
    def find_precise_x_y_min_max(self, pos_indexes, local_fits, ts, is_min):
        x_points = []
        y_points = []
        for i in range(len(pos_indexes)):
            pos_index = pos_indexes[i]
            f = local_fits[i]
            x,y = self.turn_point(ts, f, pos_index, is_min)
            x_points.append(x)
            y_points.append(y)
        return x_points, y_points
    
    # for finding minimum and maximum of first harmonic - time spectrum
    def analyse_first_harmonic_time(self, ts, fhs, is_abs, is_semi_log,
                                    is_plot_graphs, is_print_data_point_warning):
        # sort out functions for log/lin
        func = None
        if is_semi_log:
            plt.yscale('log')
            func = self.exponential_linear_func
        else:
            func = self.exponential_func
            
        max_pos_indexes, min_pos_indexes = self.coarsley_find_min_max_indices(fhs)
        
        # rough min/max pos indexes have been found
        #     now find a local fit for these
        max_pos_indexes_local_fits = self.curve_fit_min_max_regions(max_pos_indexes, ts, fhs)
        min_pos_indexes_local_fits = []
        if not is_abs:
            # need to find minima too
            min_pos_indexes_local_fits = self.curve_fit_min_max_regions(min_pos_indexes, ts, fhs)
            
        # from the local fit, find the precise location of the maxima/minima
        x_maxima, y_maxima = self.find_precise_x_y_min_max(max_pos_indexes,
                                                           max_pos_indexes_local_fits, ts, False)
        x_minima, y_minima = [],[]
        if not is_abs:
            # need to find minima too
            x_minima, y_minima = self.find_precise_x_y_min_max(min_pos_indexes,
                                                           min_pos_indexes_local_fits, ts, True)
        
        minima_fits = []
        minima_average_noise = None
        suit_xs = None
        maxima_fits, maxima_average_noise, suit_xs = self.fit(
                is_semi_log, func, True, x_maxima, y_maxima, False, is_print_data_point_warning)#is_plot_graphs, is_print_data_point_warning)
        suit_xs = array(suit_xs)
        if not is_abs:
            minima_fits, minima_average_noise, suit_min_xs = self.fit(
                    is_semi_log, func, False, x_minima, y_minima, False, is_print_data_point_warning)#is_plot_graphs, is_print_data_point_warning)
            suit_min_xs = array(suit_min_xs)
            suit_xs = concatenate([suit_xs, suit_min_xs])
            suit_xs = sort(suit_xs)
        
        noise_floor = self.find_noise_floor(is_abs, is_semi_log,
                              minima_fits, minima_average_noise,
                              maxima_fits, maxima_average_noise)
   
        damping_rate = self.find_damping_rate(is_abs, is_semi_log,
                              minima_fits, minima_average_noise,
                              maxima_fits, maxima_average_noise)
        frequency, frequency_err = self.find_average_frequency(suit_xs)
                
        if is_plot_graphs:
            self.plot_FHA_vs_time(ts, fhs)
            self.plot_suitable_point_fit(func, ts, maxima_fits, 'Maxima')
            self.plot_point_pos(x_maxima, y_maxima, True)
            self.plot_curve_fit_point_region(max_pos_indexes_local_fits, ts, max_pos_indexes, True)
            if not is_abs:
                self.plot_suitable_point_fit(func, ts, minima_fits, 'Minima')
                self.plot_point_pos(x_minima, y_minima, False)
                self.plot_curve_fit_point_region(min_pos_indexes_local_fits, ts, min_pos_indexes, False)
            plt.legend()
            plt.ioff() # This so that the windows stay open
            plt.show()
        
        return noise_floor, damping_rate, frequency, frequency_err
    
    def plot_point_pos(self, x, y, is_max):
        if is_max:
            marker, col, name = 'o', 'forestgreen', 'Maxima'
        else:
            marker, col, name = 'x', 'deepskyblue', 'Minima'
        plt.plot(x, y, marker, color=col, label=name)
    
    def plot_suitable_point_fit(self, func, ts, fits, type_name):
        # only plot if there are points to fit
        if fits is not None:
            y_fits = self.fill_array(func,fits, ts)
            plt.plot(ts, y_fits, label='Fit of Suitable {} Points'.format(type_name))
    
    def plot_curve_fit_point_region(self,local_fits, ts, pos_indexes, is_max):
        # local fits
        if is_max:
            col, type_name = 'forestgreen', 'Peak'
        else:
            col, type_name = 'deepskyblue', 'Trough'
        for i in range(len(local_fits)):
            xnew = self.turn_point_x_range(ts, pos_indexes[i])
            f = local_fits[i]
            ynew = f(xnew)
            plt.plot(xnew, ynew, linestyle='dotted', color=col)
        plt.plot(xnew, ynew, linestyle='dotted', color=col, label='Curve Fit in {} Region'.format(type_name))
        
    def plot_FHA_vs_time(self,ts,fhs):
        plt.plot(ts, fhs, color='slategray', label='Underlying Data')
        plt.xlabel("Time [Normalised]")
        plt.ylabel("First harmonic amplitude [Normalised]")
    
    def find_average_frequency(self, suit_xs):
        angular_frequencies = []
        for i in range(len(suit_xs)-1):
            time_diff = suit_xs[i+1] - suit_xs[i]
            full_time_period = 2*time_diff # with abs we see both
            frequency = 1/full_time_period
            angular_frequencies.append(2*pi*frequency) # convert to angular frequency
        # find average and standard deviation
        np_angular_frequencies = array(angular_frequencies)
        ave_frequency = mean(np_angular_frequencies)
        err = std(np_angular_frequencies)
        return ave_frequency, err
            
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