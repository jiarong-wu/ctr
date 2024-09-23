# import matplotlib as mpl


##### Make axis multiples of pi #####
# ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(
#    lambda val,pos: '{:.0g}$\pi$'.format(val/np.pi) if val !=0 else '0'
# ))
# ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base=0.5*np.pi))

from scipy.signal import hilbert
import numpy as np

############ Compute phase using Hilbert Transform ##############
def compute_phase(eta_1D):   
    analytic_signal = hilbert(eta_1D)
    phase = np.angle(analytic_signal)
    return phase