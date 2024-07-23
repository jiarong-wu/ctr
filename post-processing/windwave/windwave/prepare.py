""" Functions related to pickle. """
import sys
import os
import numpy as np
sys.path.append('/home/jiarong/research/postprocessing/functions/')
sys.path.append('/projects/DEIKE/jiarongw/jiarongw-postprocessing/jupyter_notebook/functions/')
sys.path.append('/projects/DEIKE/jiarongw/jiarongw-postprocessing/jupyter_notebook/project_specific/windwave/')
from windwave.fio import ensemble_pickle
from tqdm import tqdm

""" Helper functions """
import pickle
def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
def load_object(filename):
    with open(filename, 'rb') as input:  # Overwrites any existing file.
        obj = pickle.load(input)
    return obj

""" Helper functions (readin): 
    For the case objects, with already specified self.field_t, read in fields and store as 
    spanwise averaged time sequence of 2D array. If the pickle is available, read in the pickle;
    if not, store the pickle. (NOTICE: need to have a sub-directory of pickle_tiger/ first!)
    Output: case.ux_2D, case.uy_2D, case.f_2D (dimension time*x*z)                   
"""
def read_fields (case):
    case.ux_2D = []
    case.uy_2D = []
    case.f_2D = []

    for i in tqdm(range(0,np.size(case.field_t))):

        NSLICE = 256    
        NGRID = 512
        ux_3D = {'name':'ux', 'value':[]} # axis0 in z, axis1 in x, axis2 in y  (in the code)
        uy_3D = {'name':'uy', 'value':[]}
        f_3D = {'name':'f', 'value':[]}
        tsimu = case.field_t[i] + case.tstart
        print(tsimu)
        phasei = np.where(np.isclose(np.array(case.phase['t']), case.field_t[i]))[0][0]
        idx = case.phase['idx'][phasei]

        # Read in the fields either from pickle or from slice data
        for field in (ux_3D,uy_3D,f_3D):         
            """NOTICE: to accomodate different pickle versions"""
            picklename = case.path + 'field/' + 'pickle_tiger/' + field['name']+'_t%g' % tsimu +'.pkl'
    #             picklename = working_dir + 'field/' + 'pickle_desktop/' + field['name']+'_t%g' % t +'.pkl'
            exists = os.path.exists(picklename)
            # If the pickle is there read in the pickles
            if exists:
                field['value'] = load_object(picklename)
                print('pickle restored!')
            # If no pickle read in from the slice files and pickle dump
            if not exists:
                for sn in range (0, NSLICE-1):
                    filename = case.path + 'field/'+field['name']+'_t%g_slice%g' % (tsimu,sn)
                    snapshot = np.loadtxt(filename, dtype = np.str, delimiter='\t')
                    snapshot.reshape([NGRID,NGRID+1])
                    field['value'].append(snapshot[:,0:NGRID].astype(np.float))
                field['value'] = np.array(field['value'])
                save_object(field['value'], picklename)

            # Shift the values along x axis
            field['value'] = np.roll(field['value'], -idx, axis=1)

        case.ux_2D.append(np.average(ux_3D['value'], axis=0))
        case.uy_2D.append(np.average(uy_3D['value'], axis=0))
        case.f_2D.append(np.average(f_3D['value'], axis=0))
        
""" Helper function (readin): 
    With already specified p['t'], read in aerodynamic pressure and store as p_2D spanwise averaged time sequence of 2D array.
    The pressure here has not been processed (meaning remove the average because of the incompressible flow absolute value of pressure
    does not have meaning.) It still needs to be passed in processing_energy1 to clean the data.
    Attributes: 
        p['t']
        p['p_2D']: instantaneous pressure in x-z plane, already averaged spanwise.
        p['phat'], p['amp']: phase and amplitude of the pressure signal
"""

def read_p (case):
    case.p_2D = []
    NSLICE = 256
    NGRID = 512

    for i in tqdm(range(0, np.size(case.p['t']))):    
        pair_3D = {'name':'pair', 'value':[]}
        f_3D = {'name':'f', 'value': []}
        tsimu = case.p['t'][i] + case.tstart
        print(tsimu)
        phasei = np.where(np.isclose(np.array(case.phase['t']), case.p['t'][i]))[0][0] # Because t is a float
        idx = case.phase['idx'][phasei]

        # Read in the fields either from pickle or from slice data
        field = pair_3D
        """NOTICE: to accomodate different pickle versions"""
        picklename = case.path + 'field/' + 'pickle_tiger/' + field['name']+'_t%g' % tsimu +'.pkl'
    #             picklename = working_dir + 'field/' + 'pickle_desktop/' + field['name']+'_t%g' % t +'.pkl'
        exists = os.path.exists(picklename)
        # If the pickle is there read in the pickles
        if exists:
            field['value'] = load_object(picklename)
            print('pickle restored!')
        # If no pickle read in from the slice files and pickle dump
        if not exists:
            for sn in range (0, NSLICE-1):
                filename = case.path + 'field/'+field['name']+'_run'+'_t%g_slice%g' % (tsimu,sn)
                snapshot = np.loadtxt(filename, dtype = np.str, delimiter='\t')
                snapshot.reshape([NGRID,NGRID+1])
                field['value'].append(snapshot[:,0:NGRID].astype(np.float))
            field['value'] = np.array(field['value'])
            save_object(field['value'], picklename) 
        field['value'] = np.roll(field['value'], -idx, axis=1)

        case.p['p_2D'].append(np.average(pair_3D['value'], axis=0))


""" A go-to filter function, by default N=512, CUT=4 """
# from scipy.signal import butter,filtfilt
# def butter_lowpass_filter(data, CUT=4, N=512):
#     T = 4     # Sample Period
#     fs = 512/T      # Sample rate, Hz
#     cutoff = CUT    # desired cutoff frequency of the filter, Hz, slightly higher than actual 1.2 Hz
#     nyq = 0.5 * fs  # Nyquist Frequency
#     order = 2       # sin wave can be approx represented as quadratic
#     normal_cutoff = cutoff / nyq
#     # Get the filter coefficients 
#     b, a = butter(order, normal_cutoff, btype='low', analog=False)
#     y = filtfilt(b, a, data)
#     return y

# os.chdir('/home/jiarong/research/projects/turbulence/preliminary_cluster/stopforcing_restore_second')

# def main():
#     clock = np.arange(900, 950)
#     ensemble_pickle(clock, picklename="ensemble")
#     print(os.getcwd())

# if __name__ == "__main__":
#     main()

# print("DONE")





