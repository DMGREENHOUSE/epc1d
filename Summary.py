# -*- coding: utf-8 -*-
"""
Summary Class
"""

# Need an FFT routine, either from SciPy or NumPy
try:
    from scipy.fftpack import fft
except:
    # No SciPy FFT routine. Import NumPy routine instead
    from numpy.fft import fft

from calc_density import calc_density

class Summary:
    def __init__(self, is_print_harmonics):
        self.t = []
        self.firstharmonic = []
        self.firstharmonic_no_abs = []
        self.is_print_harmonics = is_print_harmonics
        
    def __call__(self, pos, vel, ncells, L, t):
        # Calculate the charge density
        d = calc_density(pos, ncells, L)
        
        # Amplitude of the first harmonic
        fh = 2.*abs(fft(d)[1]) / float(ncells)
        fh_no_abs = 2.*fft(d)[1] / float(ncells)
        
        if self.is_print_harmonics:
            print("Time:", t, "First:", fh)
        
        self.t.append(t)
        self.firstharmonic.append(fh)
        self.firstharmonic_no_abs.append(fh_no_abs)