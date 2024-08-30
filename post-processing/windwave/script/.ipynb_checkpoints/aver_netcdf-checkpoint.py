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
print(times)
start_time = 45

ds = read_fields_2D(path, times, NSLICE=512, NGRID=512, varlist=['ux','uy','uz','f'])
filename = path + 'span_aver.nc'
encoding = {}
for var_name in ds.data_vars:
    encoding[var_name] = {'dtype': 'float32', 'zlib': True}
ds.to_netcdf(filename, encoding=encoding)
