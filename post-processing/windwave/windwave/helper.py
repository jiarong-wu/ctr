import numpy as np
from matplotlib import pyplot as plt
import sys
sys.path.append('../../functions/')
from windwave.fio import readin
from Amplitude import Amplitude

def from_matrix(pfile):
    data = np.fromfile(pfile, dtype=np.float32)
    data = data.reshape((513,513))
    data = data[1:,1:] 
    return data

def fields(common_path, time):
    pfile = common_path + '/matrix/u%.7g.dat' %time
    u = from_matrix(pfile)
    pfile = common_path + '/matrix/f%.7g.dat' %time
    f = from_matrix(pfile)
    pfile = common_path + '/matrix/omega%.7g.dat' %time
    omega = from_matrix(pfile)
    omega_air = omega*(1-f)
    omega_water = omega*f
    u_air = u*(1-f)
    u_water = u*f
    return u_air, u_water, omega_air, omega_water

def fields_original(common_path, time):
    pfile = common_path + '/matrix/u%.7g.dat' %time
    u = from_matrix(pfile)
    pfile = common_path + '/matrix/f%.7g.dat' %time
    f = from_matrix(pfile)
    pfile = common_path + '/matrix/omega%.7g.dat' %time
    omega = from_matrix(pfile)
    return u,f,omega

def draw_field(common_path, L0, field, time, ax, absmax=800):
    image = np.rot90(field)
#     pcontour = ax.imshow(image, extent=(-L0/2,L0/2,-L0/2,L0/2), cmap='RdBu', vmax=absmax, vmin=-absmax)
    pcontour = ax.imshow(image, extent=(-0.5,0.5,-0.5,0.5), cmap='RdBu', vmax=absmax, vmin=-absmax)
    etafile = common_path + '/field/eta%.7g' %time
    eta, exists = readin(etafile, table_delimiter = ',')
    if exists:
        eta.rename(columns={'pos':'eta'}, inplace=True)
        ampl = Amplitude(eta[['x', 'eta', 'f']], 512, L0)    
    pvof = ax.plot(ampl.x_interp/L0, ampl.eta_interp/L0, color='k', linewidth=0.5)
    ax.set_title('t=%gT' %time)
    return pcontour,pvof

def interface(common_path, Npoint=512, L0=1, time=0):
    etafile = common_path + '/field/eta%.7g' %time
    eta, exists = readin(etafile, table_delimiter = ',')
    if exists:
        eta.rename(columns={'pos':'eta'}, inplace=True)
        ampl = Amplitude(eta[['x', 'eta', 'f']], Npoint, L0)   
    return ampl


'''
###############################################################################
# A class for wave properties
# Requirements:
# from scipy.optimize import fsolve
###############################################################################

'''
from scipy.optimize import fsolve
class RealWave:
    '''
    Class for calculating a set of physical properties and non-dimensional numbers
    of a monochromic wave given some quantities. All the quantities 
    are in SI units.
    '''
    
    def __init__(self, g = 9.8, sigma = 0.074, rho = 1000, rho_air = 1.225, 
                 mu = 8.9e-4, mu_air = 17.4e-6):
        '''
        Parameters
        ----------

        g : gravity acceleration (m/s^2)  
        sigma : surface tension of water/air interface (N/m)
        rho : density of water (kg/m^3)
        rho_air : density of air
        
        self.k : wave number (1/m)
        self.omega : wave frequency (1/s)
        self.c : phase speed (m/s)
        self.wl : wavelength (m)
        self.Bo : Bond number 
        
        '''
        self.g, self.sigma, self.rho, self.rho_air, self.mu, self.mu_air = \
        g, sigma, rho, rho_air, mu, mu_air
        self.k, self.omega, self.c, self.wl, self.Bo, self.Re_wave, self.Re_air = 0, 0, 0, 0, 0, 0, 0
        
    def k2omega(self,k):
        self.k = k
        # Gravity-capillary wave dispersion relation
        self.omega = (self.g*self.k + self.sigma*self.k**3/self.rho)**0.5
        self.c = self.omega/self.k
        self.wl = 2*np.pi/self.k
        self.Bo =  (self.rho-self.rho_air)*self.g/self.sigma/self.k**2
#         print("Given k = %g (1/m), calculated omega = %g (1/s), period = %g (s), phase speed c = %g (m/s), wavelength = %g (m), Bo = %g" 
#               %(self.k, self.omega, 2*np.pi/self.omega, self.c, self.wl, self.Bo))

    # Implicit function of w(k)
    def omega2k(self,omega):
        self.omega = omega
        k = fsolve(lambda k : (self.g*k + self.sigma*k**3/self.rho)**0.5 - omega, 0)
        self.k = k[0]
        self.c = self.omega/self.k
        self.wl = 2*np.pi/self.k
        self.Bo =  (self.rho-self.rho_air)*self.g/self.sigma/self.k**2
#         print("Given omega = %g (1/s), calculated k = %g (1/m), phase speed c = %g (m/s), wavelength = %g (m), Bo = %g" 
#               %(self.omega, self.k, self.c, self.wl, self.Bo))
              
    # If Bond number is given instead of k
    def Bo2k(self,Bo):
        self.Bo = Bo
        self.k = ((self.rho-self.rho_air)*self.g/Bo/self.sigma)**0.5
        self.wl = 2*np.pi/self.k
        self.omega = (self.g*self.k + self.sigma*self.k**3/self.rho)**0.5
        self.c = self.omega/self.k
#         print("Given Bo = %g, calculated lambda = %g (m), k = %g (1/m), omega = %g (1/s), phase speed c = %g (m/s)" 
#               %(self.Bo, self.wl, self.k, self.omega, self.c))
     
    def Re(self, UstarRatio = 0, L0 = 1.):
        self.Re_wave = self.rho*self.c*(2*np.pi/self.k)/self.mu
        self.Re_air = self.rho_air*self.c*UstarRatio*(L0/2.)/self.mu_air
        c_simu = (1/2/np.pi*(1+1/self.Bo))**0.5
        self.Re_nominal = self.Re_wave/c_simu
#         print("Re_wave = %g, Re_air = %g, Re in older version = %g" %(self.Re_wave, self.Re_air, self.Re_nominal))