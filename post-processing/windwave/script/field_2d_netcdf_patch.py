"""
Script Name: field_2d_netcdf_patch.py

Description:
    This script patches missing time snapshots to netcdf file assembled using field_2d_netcdf.py
    filling the NaNs with new *_2d_avg_t*.bin files and
    save to new file named span_aver_patched.nc

Usage:
    Manually give the following:
        path: path of the case to patch
        path_patch: where the patching field is stored and time instances to patch
    Run it with "python field_2d_netcdf_patch.py" in command line.

Dependence:
    windwave.fio where the read_fields_2D function is defined
    xarray and numpy     
"""

import xarray as xr
from windwave.fio import read_fields_2D

path = '/home/ctrsp-2024/jiarongw/outputs/CU8/'
ds = xr.load_dataset(f'{path}span_aver.nc')

path_patch = '/home/ctrsp-2024/jiarongw/outputs/CU8/missing/'
ds_patch = read_fields_2D(path_patch, [52.2], NSLICE=512, NGRID=512, varlist=['ux','uy','uz','f'])

ds = ds.where(~np.isnan(ds.ux), ds_patch.ux.values)
ds = ds.where(~np.isnan(ds.uy), ds_patch.uy.values)
ds = ds.where(~np.isnan(ds.uz), ds_patch.uz.values)
ds = ds.where(~np.isnan(ds.f), ds_patch.f.values)

filename = path + 'span_aver_patched.nc'
encoding = {}
for var_name in ds.data_vars:
    encoding[var_name] = {'dtype': 'float32', 'zlib': True}
ds.to_netcdf(filename, encoding=encoding)
