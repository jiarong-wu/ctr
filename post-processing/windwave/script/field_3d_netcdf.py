"""
Script Name: field_3d_netcdf.py

Description:
    This script read Basilisk raw output of slices (either along Z or Y) to construct 
    the full 3D field. Files are read from specified case _path_. 
    Fields given by _varlist_ are aggreated into one xarray dataset,
    Raw 2D slices should be named and stored as '${path}/field/*_t*_zslice*' or '${path}/field/*_t*_yslice*'.
    Otherwise modify the windwave.fio accordingly.
    
Usage:
    Manually give the following:
        path: path of the case
        time: an array of available time snapshots
        start_time: the time when wave starts moving
    Run it with "python field_2d_netcdf.py" in command line.

Dependence:
    windwave.fio where the read_fields_2D function is defined
    xarray and numpy     
"""
import argparse
import glob
import os
import numpy as np
import xarray as xr
import gc
from scipy.interpolate import griddata
from windwave.fio import read_fields_3D, read_eta
from windwave.utils import compute_phase

#### Manually set the case path and some meta data, although we could read them from command line ####
DATAPATH = '/scratch/jw8736/ctr/outputs/'

# path = DATAPATH + 'NWP_ZPG_CU2/'
# times = np.arange(45,59)

# path = DATAPATH + 'NWP_ZPG_CU4/'
# times = np.arange(41,59)

path = DATAPATH + 'NWP_ZPG_CU8/'
times = np.arange(50,51)

# path = DATAPATH + 'CU4/'
# times = np.arange(56,59)

# path = DATAPATH + 'CU8/'
# times = np.arange(48,61)

# A list of meta data
start_time = 45
ak = 0.2
k = 4
L0 = 2*np.pi

if __name__ == "__main__":
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="Assemble the 3D field and save to netcdf, either from scratch or from an already assembled dataset missing eta.")
    parser.add_argument('--mode', type=int, required=True, help="From scratch (0) or from already assembled dataset missing eta (1)")
    # Parse the arguments
    args = parser.parse_args()

    print(path)
    
    #### If starting from scratch and reading in fields and eta together #####
    if args.mode == 0:
        for t in times:
            print('Processing t=%g... from scratch' %t)
        
            # STEP1: read in 3D fields
            ds = read_fields_3D(path, t=t, NSLICE=512, NGRID=512, varlist=['ux','uy','uz','f'], SLICE='Z', SAVE=False)
            
            ds.attrs['start_time'] = start_time
            ds.attrs['ak'] = ak
            ds.attrs['k'] = k
            ds.attrs['L0'] = L0
        
            # STEP2: read in eta, compute phase (spanwise averaged), and append eta and phase fields to dataset 
            xtile, ztile = np.meshgrid(ds.x.values, ds.z.values) # xy indexing meaning that z rows and x columns
            
            if t > ds.attrs['start_time']:
                # read file of columns of x, z, eta, and slope
                etaname = path + 'eta/eta_t%g' %t
            elif t <= ds.attrs['start_time']:
                # initial slope of stationary eta
                etaname = path + 'eta/eta_initial'
                
            snapshot = read_eta(etaname)    
            # interpolate onto uniform x-z grid
            eta = griddata((snapshot.x, snapshot.z), snapshot.pos, (xtile, ztile), method='linear')
         
            # append eta and phase fields
            ds['eta'] = (['z','x'], eta)
            eta_1D = ds['eta'].mean(dim='z') - ds['eta'].mean(dim=['z','x'])
            ds['phase'] = (('x'), compute_phase(eta_1D))
        
            # Writing to new file 
            encoding = {}
            for var_name in ds.data_vars:
                encoding[var_name] = {'dtype': 'float32', 'zlib': True}
        
            file = path + 'netcdf/field_eta_t%g.nc' %t
            ds.to_netcdf(file, encoding=encoding)

            del ds # this is needed for memory 
            gc.collect()
    
    
    #### If file has already been saved to netcdf/field_t*.nc using read_fields_3D with SAVE=True #####
    #### Additional read-in of eta and compute phase #####
    if args.mode == 1:
        
        files = sorted(glob.glob(os.path.join(path, 'netcdf/field_t*.nc')))  # Use wildcard to find all matching files
        
        for file in files:
            t = int(os.path.basename(file).split('_t')[1].split('.nc')[0])
            ds = xr.open_dataset(file)
            
            ds.attrs['start_time'] = start_time
            ds.attrs['ak'] = ak
            ds.attrs['k'] = k
            ds.attrs['L0'] = L0
            
            print('Processing t=%g... only appending eta' %t)
            xtile, ztile = np.meshgrid(ds.x.values, ds.z.values) # xy indexing meaning that z rows and x columns
            
            if t > ds.attrs['start_time']:
                # read file of columns of x, z, eta, and slope
                etaname = path + 'eta/eta_t%g' %t
            elif t <= ds.attrs['start_time']:
                # initial slope of stationary eta
                etaname = path + 'eta/eta_initial'
                
            snapshot = read_eta(etaname)    
            # interpolate onto uniform x-z grid
            eta = griddata((snapshot.x, snapshot.z), snapshot.pos, (xtile, ztile), method='linear')
        
            # append eta and phase fields
            ds['eta'] = (['z','x'], eta)
            eta_1D = ds['eta'].mean(dim='z') - ds['eta'].mean(dim=['z','x'])
            ds['phase'] = (('x'), compute_phase(eta_1D))
        
            # Writing to new file 
            encoding = {}
            for var_name in ds.data_vars:
                encoding[var_name] = {'dtype': 'float32', 'zlib': True}
        
            newfile = file.replace("field_", "field_eta_")
            ds.to_netcdf(newfile, encoding=encoding)
    
    print('Done!')


#### If we want to delete old files ####
# try:
#     os.remove(file_path)
#     print(f"File {file_path} deleted successfully.")
# except FileNotFoundError:
#     print(f"File {file_path} not found.")
# except PermissionError:
#     print(f"Permission denied: unable to delete {file_path}.")
# except Exception as e:
#     print(f"Error occurred while deleting the file: {e}")