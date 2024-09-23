"""
Script Name: field_2d_netcdf.py

Description:
    This script read Basilisk raw output of span-wise averaged 2D fields
    from specified case _path_, aggregate all fields given by _varlist_ into one xarray dataset,
    concatenate along time axis, and save the concatenated xarray dataset to the same path
    with the file name "span_aver.nc", with some compression.
    Raw 2D fields should be named as f_2d_avg_t45.bin, ux_2d_avg_t45.bin, ... and stored in path/field/
    
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



from windwave.fio import read_fields_2D
import xarray as xr
import numpy as np

# path = '/home/ctrsp-2024/jiarongw/outputs/CU2/'
# times = np.arange(40,52,0.1)
# start_time = 40

# path = '/home/ctrsp-2024/jiarongw/outputs/CU4/'
# times = np.arange(38.5,59,0.1)
# start_time = 40

# path = '/home/ctrsp-2024/jiarongw/outputs/CU8/'
# times = np.arange(40,60.6,0.1)
# start_time = 40

# path = '/home/ctrsp-2024/jiarongw/outputs/NWP_ZPG_CU2/'
# times = np.arange(45,58,0.1)
# start_time = 45

# path = '/home/ctrsp-2024/jiarongw/outputs/NWP_ZPG_CU4/'
# times = np.arange(41.5,58,0.1)
# start_time = 45

# path = '/home/ctrsp-2024/jiarongw/outputs/NWP_ZPG_CU8/'
# times = np.arange(45,59,0.1)
# start_time = 45

# path = '/home/ctrsp-2024/jiarongw/outputs/NWP_ZPG_CU8_ensem2/'
# times = np.arange(45,59,0.1)
# start_time = 45

# path = '/home/ctrsp-2024/jiarongw/outputs/NWP_ZPG_CU8_ensem3/'
# times = np.arange(45,58,0.1)
# start_time = 45

path = '/home/ctrsp-2024/jiarongw/outputs/test_NWP_precursor/'
times = np.arange(45,46.6,0.1)
start_time = 45

print(times)

ds = read_fields_2D(path, times, NGRID=512, varlist=['ux','uy','uz','f'])
filename = path + 'span_aver.nc'
encoding = {}
for var_name in ds.data_vars:
    encoding[var_name] = {'dtype': 'float32', 'zlib': True}
ds.to_netcdf(filename, encoding=encoding)
