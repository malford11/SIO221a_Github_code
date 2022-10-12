#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import xarray as xr
import glob

def GetPierData(variables):
    """
    This function creates a dataset of all requested variables at the Scripps Pier station from 2005-2021.
    
    variables should be in string format, i.e. 'temperature', 'pressure', etc. 
    multiple variables can be specified as ['temperature','pressure', ...]
    
    """

    filenames = glob.glob('scripps_pier*') # Collecting all filenames

    ds_list = [] # Empty list to put each year's dataset in
    for f in filenames:
        # Looping over all files
        ds = xr.open_dataset(f) # xarray package reads netCDF files
        ds = ds[variables] # selecting only variables specified in function input
        time_fixed = [np.datetime64(ti) for ti in ds.time.data] # converting time coordinate to datetime64 format
        ds['time'] = time_fixed
        ds_list.append(ds) # Appending each year to ds_list
        
    ds_combined = xr.concat(ds_list, dim='time') # Combining all years into one big dataset
    
    return ds_combined

