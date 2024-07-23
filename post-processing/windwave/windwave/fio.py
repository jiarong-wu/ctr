import os
import numpy as np
import xarray as xr


############ Read in 2D span-wise averaged data, along time axis ##############

def read_fields_2D (path, times, NSLICE=256, NGRID=512, varlist=['ux','uy','uz','f']):
    
    filename = path + 'netcdf/span_aver.nc' 

    if os.path.exists(filename):
        print('NetCDF file exist!')
        ds = xr.open_dataset('filename')
    
    else:
        # if not construct new xarray dataset
        x = np.linspace(-np.pi, np.pi, NGRID, endpoint=False) + 2*np.pi/NGRID/2
        y = np.linspace(0, 2*np.pi, NGRID, endpoint=False) + 2*np.pi/NGRID/2
        
        # TODO: write read in eta file
        # ...

        # read in 2D span-wise averaged field
        for i, var in enumerate(varlist):   
            field = []
            for t in times:
                name = path + 'field/' + var + '_2d_avg_t%g.bin' % (t)
                avg = np.fromfile(name, dtype=np.float64)
                avg = avg.reshape([NGRID,NGRID])
                field.append(avg)
            field = np.array(field)
            if i == 0:
                ds = xr.Dataset({var: (['t','x','y'], field)}, coords={'t':times, 'x':x, 'y':y})
            else:
                ds = ds.assign(**{var: (['t','x','y'], field)})

        # Writing to file with some compression
        # encoding = {}
        # for var_name in ds.data_vars:
        #     encoding[var_name] = {'dtype': 'float32', 'zlib': True}
        # ds.to_netcdf(filename, encoding=encoding)

        # Shift the values along x axis
        # field['value'] = np.roll(field['value'], -idx, axis=1)
        return ds

############ Old code ###############

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


################## old code #####################
import os
import pandas as pd
import numpy as np
import pickle
from tqdm import tqdm

def readin(filename, table_delimiter = ',', table_headers = None, skipn = None):
    '''
    A function used to read in field contained data file in the form of 
    dataframe. The warning of specifying dtype is to be ignored because of
    the current way of first reading in data then convert to numerics.
    Some helpful features are:
    1. Check if the file is damaged
    2. avoid corrupted lines by forcing converting to numeric and drop NAN
    
    filename: string
        The name of the data file. (Not supposed to be iterative. Iterating is performed by a 
        higher level function.)
    table_delimiter: string, optional
        Either '' or ',', default ','
    table_headers: list of keys, optional
        Names of the field attributes. Default is None.    
    skipn : int, optional
        Skip every skipn rows. Apply when the input data is too dense.
    '''
    
    # An implementation of droppping lines according to the return value of logic
#     def logic(index):
#         if index % 100 == 0:
#             return False
#         return True
#         energy = pd.read_table(filename, delimiter = ' ', skiprows= lambda x: logic(x), error_bad_lines=False)
    # An implementation of dropping duplicate in a certain column
#         data = data.drop_duplicates(subset=['t'], keep='last')

    exists = os.path.exists(filename)
    if not exists:
        print(filename + ' cannot be read!')        
    if exists:
        # Skip rows or not
        if skipn!=None:   
            def logic(index):
                if index % skipn == 0:
                   return False
                return True
            data = pd.read_table(filename, delimiter = table_delimiter, skiprows= lambda x: logic(x), error_bad_lines=False)
        else:
            data = pd.read_table(filename, delimiter = table_delimiter, names = table_headers, error_bad_lines=False)
        # Convert to numerical
        columns = list(data) 
        for i in columns: 
            data[i] = pd.to_numeric(data[i],errors='coerce')
        data = data.dropna()
        data = data.reset_index(drop=True)
    else:
        data = None
    return data, exists

def put_together(filename, core_number, table_delimiter = ',', table_headers = None):
    '''
    A function used to put data files from different cores together.
    
    filename: string
        The common part of the name of the data file. (Not supposed to be iterative.)
    core_number: int
        Number of cores used during the parallel computing.
    table_delimiter, table_headers: same as above
      
    '''
    filename_number = filename + '_0.dat' 
    data = pd.read_table(filename_number, names = table_headers, delimiter = table_delimiter)
    for i in range (0, core_number):
        filename_number = filename + '%g.dat' % i
        data_append = pd.read_table(filename_number, names = table_headers, delimiter = table_delimiter)
        data = data.append(data_append, ignore_index=True)
    return data

def ensemble_pickle(clock, headers=['x', 'y', 'z', 'u.x', 'u.y', 'u.z', 'omega'], 
                    filename_common = './field_direct', picklename = 'ensemble', sort_and_halve = True):
    '''
    Read data across different dump files. Better to be called from a separate script instead
    of inside a notebook cell. See /postprocessing/turbulence/prepare.py
    Put together in list field; record respective time in list tseries;
    count number of points in list ptnumber(for checking if maximum level is reached).
    
    clock: list of time
    headers: keys of data file headers
    filename_common: string of the file name without time postsuffix
    picklename: string of the output pickle file name
    sort_and_halve: boolean to perform extra data clean, optional, default True    
    '''
    
    field = []
    tseries = []
    ptnum = []
    for t in tqdm(clock):
        filename = filename_common + '%g' % t
        snapshot = readin(filename, table_headers = headers)
        # Sort dataframe by y value, only taking the upper half, optional
        if sort_and_halve:
            snapshot = snapshot.sort_values(by = ['x','y','z'])
            snapshot = snapshot.loc[snapshot.y >= 0]
            # Reset the index which is important! AlWAYS reset the fucking index!!!
            snapshot = snapshot.reset_index(drop=True)
        field.append(snapshot)
        tseries.append(t)
        ptnum.append(len(snapshot))
    # After reading all the data, pickle them  
    outfile = open(picklename,'wb')
    data = (field, tseries, ptnum)
    pickle.dump(data,outfile)