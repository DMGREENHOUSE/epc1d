# -*- coding: utf-8 -*-
"""
Class for a Data Run
"""
from Summary import Summary
from Plot import Plot
from numpy import arange, concatenate, linspace, floor, array, pi
from numpy import sin, cos, sqrt, random



# Need an FFT routine, either from SciPy or NumPy
try:
    from scipy.fftpack import fft, ifft
except:
    # No SciPy FFT routine. Import NumPy routine instead
    from numpy.fft import fft, ifft

from calc_density import calc_density

class Data_Run:
    def __init__(self, is_Landau=True, L=None, ncells=None, npart=None,
                 time_length = 20, num_time_steps = 50,
                 is_plot_animation = False, is_print_harmonics = True,
                 iteration_num=0):
        # Generate initial condition
        #
        self.is_Landau = is_Landau
        if is_Landau:
            # Landau damping
            self.L = 4.*pi if L == None else L
            self.ncells = 20 if ncells == None else ncells
            self.npart = 1000 if npart == None else npart
            
            pos, vel = self.landau(self.npart, self.L)
        else:
            # 2-stream instability
            self.L = 100 if L == None else L
            self.ncells = 20 if ncells == None else ncells
            self.npart = 10000 if npart == None else npart
            pos, vel = self.twostream(self.npart, self.L, 3.)
        
        self.time_length = time_length
        self.num_time_steps = num_time_steps
        self.iteration_num = iteration_num
        # Create some output classes
        p = Plot(pos, vel, self.ncells, self.L, is_plot_animation) # This displays an animated figure
        self.s = Summary(is_print_harmonics)                 # Calculates, stores and prints summary info
        
        # Run the simulation
        pos, vel = self.run(pos, vel, self.L, self.ncells, 
                       out=[p, self.s],                      # These are called each output
                       output_times=linspace(0.,self.time_length,self.num_time_steps)) # The times to output
        
    def get_summary(self):
        return self.s
        
    def get_type(self):
        this_type = None
        if self.is_Landau:
            this_type = 'LD'
        else:
            this_type = '2SI'
        return this_type
    def get_iteration_num(self):
        return self.iteration_num
    def get_L(self):
        return self.L
        
    def get_ncells(self):
        return self.ncells
        
    def get_npart(self):
        return self.npart
        
    def get_time_length(self):
        return self.time_length
        
    def get_num_time_steps(self):
        return self.num_time_steps
    
    def run(self, pos, vel, L, ncells=None, out=[], output_times=linspace(0,20,100), cfl=0.5):
        
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
                f = self.rk4step(self.pic, f, dt, args=(ncells, L))
                time += dt
                
            # Extract position and velocities
            pos = ((f[0:nparticles] % L) + L) % L
            vel = f[nparticles:]
            
            # Send to output functions
            for func in out:
                func(pos, vel, ncells, L, time)
            
        return pos, vel
    
    def periodic_interp(self, y, x):
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
    
    def fft_integrate(self, y):
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
       
    
    def pic(self, f, ncells, L):
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
        E = -self.fft_integrate(rho)*dx
        
        # Interpolate E field at particle locations
        accel = -self.periodic_interp(E, pos/dx)
    
        # Put back into a single array
        return concatenate( (vel, accel) )

    def rk4step(self, f, y0, dt, args=()):
        """ Takes a single step using RK4 method """
        k1 = f(y0, *args)
        k2 = f(y0 + 0.5*dt*k1, *args)
        k3 = f(y0 + 0.5*dt*k2, *args)
        k4 = f(y0 + dt*k3, *args)
    
        return y0 + (k1 + 2.*k2 + 2.*k3 + k4)*dt / 6.
    
    ####################################################################
    # 
    # Functions to create the initial conditions
    #
    
    def landau(self, npart, L, alpha=0.2):
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
    
    def twostream(self, npart, L, vbeam=2):
        # Start with a uniform distribution of positions
        pos = random.uniform(0., L, npart)
        # Normal velocity distribution
        vel = random.normal(0.0, 1.0, npart)
        
        np2 = int(npart / 2)
        vel[:np2] += vbeam  # Half the particles moving one way
        vel[np2:] -= vbeam  # and half the other
        
        return pos,vel
    
    ####################################################################
         
    ####################################################################