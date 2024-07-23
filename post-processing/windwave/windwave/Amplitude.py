# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 09:42:48 EDT 2019
@author: jiarongw
"""

import numpy as np
import pandas as pd
import scipy.interpolate
import os
import pickle
import matplotlib.pyplot as plt
from scipy import signal

'''
###############################################################################
# A class for dealing with amplitude
###############################################################################
'''

class Amplitude:
    '''
    Class for analyzing the interface elevation. Instantiation takes the eta data.
    
    self._eta_data: dataframe
        The complete dataframe, sorted by value of x.
    self.eta: 1D array
        Only the amplitude info without coordinate.
    self.stdev: float
        Standard deviation of amplitude. Calculated by std = sqrt(mean(abs(x - x.mean())**2)).
    self.phase: float
        Given a certian elevation profile, there should be a 
            
    '''
    
    def __init__(self, eta_data, N = 512, L0 = 1):
        '''
        eta_data : dataframe
            A set of data containing the position x, interface elevation eta, and scalar fraction
            field f. The corresponding keywords are ['x', 'eta', 'f']. f is used for filtering out 
            stand alone little points.
            Later a 3D feature might be added.
        '''   
        # A clean-up process
        # First filter out points with extreme values of f (tolerance can be adjusted)
        # Then sort the array according to position x
        tol = 1e-5
        self._eta_data = eta_data.loc[(eta_data.f > tol) & ((1-eta_data.f) > tol)]
        self._eta_data = self._eta_data.sort_values(by = ['x'])
        self.eta = eta_data['eta'].values
        self.x = eta_data['x'].values
        # Compute slope
#         self.slope = np.gradient(self.eta, self.x) 
        # Some polyfit ...
        if np.any(self.eta): # prevent an all zero array as it happens sometimes
            self.x_interp = np.linspace(-L0/2,L0/2,N,endpoint=False)+L0/N/2
            self.eta_interp = np.interp(self.x_interp, self._eta_data.x, self._eta_data.eta)
            self.stdev_interp = np.std(self.eta_interp)
            self.stdev = np.std(self.eta)
        else:
            self.stdev = 0        
            self.stdev_interp = 0
    
    # A plotting function for amplitude profile
    # Use non-dimensionalized amplitude ak or not?
    def plot_amplitude(self, ax, label_choice = None, color_choice = None):
        ax.plot(self._eta_data.x, self._eta_data.eta, label = label_choice, color = color_choice)
        ax.set_xlabel('x position')
        ax.set_ylabel('amplitude')
    
    #A plotting function for the spectrum
    def plot_spectrum(self, ax, cutoff = 8, label_choice = None, color_choice = None):
        # Plot the magnitude of wave spectrum
        '''
        cutoff: integer, optional
            Cut off the plotting at the nth wavenumber since higher frequency 
            contains no energy. Default is 8.
        '''
        ax.plot(self.spectrum[0][0:cutoff], abs(self.spectrum[1])[0:cutoff], 
            label = label_choice, color = color_choice)
        ax.set_xlabel('wavenumber')
        ax.set_ylabel('Y(Wavenumber)')

        
    def phase(self):
        phase = 0
        '''
        There is a uniquely determined phase for a given amplitude profile
        '''

    def FFT(self, domain = None, Fs = 2048):
        '''
        domain_width: tuple, optional
            A tuple with two component, starting point and end point. 
            If not specified, take the extrema of self._eta_data.x
        Fs: integer, optional
            Sampling rate. Power of 2. 
        
        Reference: https://www.ritchievink.com/blog/2017/04/23/understanding-the-fourier-transform-by-example/
        
        '''
        if domain is None:
            domain = (self._eta_data.x.min(), self._eta_data.x.max())
        Xs = 1./Fs; # sampling interval
        x = np.arange(domain[0], domain[1], Xs) # space vector
        y = np.interp(x, self._eta_data.x, self._eta_data.eta)
        # In frequency domain, the highest is sampling frequency Fs, the 
        # resolution is determined by sample size: \deltaf = fs/N

        N = len(y) # sample size, may not be Fs depending on domain size
        wavenumber = np.linspace(0, Fs, N) 
        # An alternative way
        # k = np.arange(n)
        # T = n/Fs
        # wavenumber = k/T 
        wavenumber = wavenumber[0:int(N/2)] # one side frequency range
        Y = np.fft.fft(y)/N # fft computing and normalization
        Y = Y[0:int(N/2)]
        self.spectrum = (wavenumber, Y)
    


