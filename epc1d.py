#!/usr/bin/env python
#
# Electrostatic PIC code in a 1D cyclic domain

from numpy import arange, concatenate, zeros, linspace, floor, array, pi
from numpy import sin, cos, sqrt, random, histogram

from scipy import interpolate
from numpy import append, argmin, argmax

import matplotlib.pyplot as plt # Matplotlib plotting library

try:
    import matplotlib.gridspec as gridspec  # For plot layout grid
    got_gridspec = True
except:
    got_gridspec = False

# Need an FFT routine, either from SciPy or NumPy
try:
    from scipy.fftpack import fft, ifft
except:
    # No SciPy FFT routine. Import NumPy routine instead
    from numpy.fft import fft, ifft

def rk4step(f, y0, dt, args=()):
    """ Takes a single step using RK4 method """
    k1 = f(y0, *args)
    k2 = f(y0 + 0.5*dt*k1, *args)
    k3 = f(y0 + 0.5*dt*k2, *args)
    k4 = f(y0 + dt*k3, *args)

    return y0 + (k1 + 2.*k2 + 2.*k3 + k4)*dt / 6.

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
    for p in position / dx:    # Loop over all the particles, converting position into a cell number
        plower = int(p)        # Cell to the left (rounding down)
        offset = p - plower    # Offset from the left
        density[plower] += 1. - offset
        density[(plower + 1) % ncells] += offset
    # nparticles now distributed amongst ncells
    density *= float(ncells) / float(nparticles)  # Make average density equal to 1
    return density

def periodic_interp(y, x):
    """
    Linear interpolation of a periodic array y at index x
    
    Input

    y - Array of values to be interpolated
    x - Index where result required. Can be an array of values
    
    Output
    
    y[x] with non-integer x
    """
    ny = len(y)
    if len(x) > 1:
        y = array(y) # Make sure it's a NumPy array for array indexing
    xl = floor(x).astype(int) # Left index
    dx = x - xl
    xl = ((xl % ny) + ny) % ny  # Ensures between 0 and ny-1 inclusive
    return y[xl]*(1. - dx) + y[(xl+1)%ny]*dx

def fft_integrate(y):
    """ Integrate a periodic function using FFTs
    """
    n = len(y) # Get the length of y
    
    f = fft(y) # Take FFT
    # Result is in standard layout with positive frequencies first then negative
    # n even: [ f(0), f(1), ... f(n/2), f(1-n/2) ... f(-1) ]
    # n odd:  [ f(0), f(1), ... f((n-1)/2), f(-(n-1)/2) ... f(-1) ]
    
    if n % 2 == 0: # If an even number of points
        k = concatenate( (arange(0, n/2+1), arange(1-n/2, 0)) )
    else:
        k = concatenate( (arange(0, (n-1)/2+1), arange( -(n-1)/2, 0)) )
    k = 2.*pi*k/n
    
    # Modify frequencies by dividing by ik
    f[1:] /= (1j * k[1:]) 
    f[0] = 0. # Set the arbitrary zero-frequency term to zero
    
    return ifft(f).real # Reverse Fourier Transform
   

def pic(f, ncells, L):
    """ f contains the position and velocity of all particles
    """
    nparticles = len(f)/2     # Two values for each particle
    nparticles_int = int(nparticles)# convert to integer
    pos = f[0:nparticles_int] # Position of each particle
    vel = f[nparticles_int:]      # Velocity of each particle

    dx = L / float(ncells)    # Cell spacing

    # Ensure that pos is between 0 and L
    pos = ((pos % L) + L) % L
    
    # Calculate number density, normalised so 1 when uniform
    density = calc_density(pos, ncells, L)
    
    # Subtract ion density to get total charge density
    rho = density - 1.
    
    # Calculate electric field
    E = -fft_integrate(rho)*dx
    
    # Interpolate E field at particle locations
    accel = -periodic_interp(E, pos/dx)

    # Put back into a single array
    return concatenate( (vel, accel) )

####################################################################

def run(pos, vel, L, ncells=None, out=[], output_times=linspace(0,20,100), cfl=0.5):
    
    if ncells == None:
        ncells = int(sqrt(len(pos))) # A sensible default

    dx = L / float(ncells)
    
    f = concatenate( (pos, vel) )   # Starting state
    nparticles = len(pos)
    
    time = 0.0
    for tnext in output_times:
        # Advance to tnext
        stepping = True
        while stepping:
            # Maximum distance a particle can move is one cell
            dt = cfl * dx / max(abs(vel))
            if time + dt >= tnext:
                # Next time will hit or exceed required output time
                stepping = False
                dt = tnext - time
            #print "Time: ", time, dt
            f = rk4step(pic, f, dt, args=(ncells, L))
            time += dt
            
        # Extract position and velocities
        pos = ((f[0:nparticles] % L) + L) % L
        vel = f[nparticles:]
        
        # Send to output functions
        for func in out:
            func(pos, vel, ncells, L, time)
        
    return pos, vel

####################################################################
# 
# Output functions and classes
#

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
        

class Summary:
    def __init__(self, is_print_harmonics):
        self.t = []
        self.firstharmonic = []
        self.is_print_harmonics = is_print_harmonics
        
    def __call__(self, pos, vel, ncells, L, t):
        # Calculate the charge density
        d = calc_density(pos, ncells, L)
        
        # Amplitude of the first harmonic
        fh = 2.*abs(fft(d)[1]) / float(ncells)
        
        if self.is_print_harmonics:
            print("Time:", t, "First:", fh)
        
        self.t.append(t)
        self.firstharmonic.append(fh)

####################################################################
# 
# Functions to create the initial conditions
#

def landau(npart, L, alpha=0.2):
    """
    Creates the initial conditions for Landau damping
    
    """
    # Start with a uniform distribution of positions
    pos = random.uniform(0., L, npart)
    pos0 = pos.copy()
    k = 2.*pi / L
    for i in range(10): # Adjust distribution using Newton iterations
        pos -= ( pos + alpha*sin(k*pos)/k - pos0 ) / ( 1. + alpha*cos(k*pos) )
        
    # Normal velocity distribution
    vel = random.normal(0.0, 1.0, npart)
    
    return pos, vel

def twostream(npart, L, vbeam=2):
    # Start with a uniform distribution of positions
    pos = random.uniform(0., L, npart)
    # Normal velocity distribution
    vel = random.normal(0.0, 1.0, npart)
    
    np2 = int(npart / 2)
    vel[:np2] += vbeam  # Half the particles moving one way
    vel[np2:] -= vbeam  # and half the other
    
    return pos,vel

####################################################################

def is_pos_diff(vals, test_index):
    is_this_pos_grad = None
    diff = vals[test_index+1]-vals[test_index]
    if diff >= 0:
        is_this_pos_grad = True
    else:
        is_this_pos_grad = False
    return is_this_pos_grad

def turn_point(xs, interp, target_index, is_min):
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
def analyse_first_harmonic_time(ts, fhs):
    max_pos_indexes = []
    min_pos_indexes = []
    is_prev_pos_grad = is_pos_diff(fhs, 0)
    for i in range(len(fhs)-1):
        # check if this one is a positive gradient
        is_this_pos_grad = is_pos_diff(fhs, i)
        # check for change from previous
        if is_prev_pos_grad and not is_this_pos_grad:
            max_pos_indexes.append(i)
        if not is_prev_pos_grad and is_this_pos_grad:
            min_pos_indexes.append(i)
        is_prev_pos_grad = is_this_pos_grad
    f = interpolate.interp1d(ts, fhs, kind='quadratic')
    # do the first one
    x,y = turn_point(ts, f, 1, False)
    plt.plot(x, y, 'o', color='k')
    for min_index in min_pos_indexes:
        x,y = turn_point(ts, f, min_index, True)
        plt.plot(x, y, 'x', color='r')
    for max_index in max_pos_indexes:
        x,y = turn_point(ts, f, max_index, False)
        plt.plot(x, y, 'o', color='k')
    
    
    xnew = arange(ts[0], ts[-1], 0.05)
    ynew = f(xnew)   # use interpolation function returned by `interp1d`
    plt.plot(ts, fhs, 'o', xnew, ynew, '-')
    plt.plot(ts, fhs)
    plt.show()
     
####################################################################

if __name__ == "__main__":
    # Generate initial condition
    #
    if False:
        # 2-stream instability
        L = 100
        ncells = 20
        pos, vel = twostream(10000, L, 3.)
    else:
        # Landau damping
        L = 4.*pi
        ncells = 20
        npart = 1000
        pos, vel = landau(npart, L)
    is_plot_animation = False
    is_print_harmonics = True
    # Create some output classes
    p = Plot(pos, vel, ncells, L, is_plot_animation) # This displays an animated figure
    s = Summary(is_print_harmonics)                 # Calculates, stores and prints summary info
    
    # Run the simulation
    pos, vel = run(pos, vel, L, ncells, 
                   out=[p, s],                      # These are called each output
                   output_times=linspace(0.,20,50)) # The times to output
    
    # Summary stores an array of the first-harmonic amplitude
    analyse_first_harmonic_time(s.t, s.firstharmonic)
    # Make a semilog plot to see exponential damping
    plt.figure()
    plt.plot(s.t, s.firstharmonic)
    plt.xlabel("Time [Normalised]")
    plt.ylabel("First harmonic amplitude [Normalised]")
    plt.yscale('log')
    
    plt.ioff() # This so that the windows stay open
    plt.show()
    
    
