def Get46047DataMY(path,year_start=None,year_end=None):
    # This function returns a multi-year wind speed, wave height, water temperature, and air temperature 
    # from  the 46047 buoy.
    # (https://www.ndbc.noaa.gov/historical_data.shtml)
    
    # The path must be a directory containing all the 46047 buoy data from 2015-2020
    # year_start is the start year of interest (must be between 2015 and 2020)
    # year_end is the end year of interest (must be between 2015 and 2020)
    # Default is year_start and year_start equal to None (all the 2015-2020 record will be loaded)
    # CDL SIO211A 10/13/2022

    # Important packages
    import datetime                   # Manage time vectors
    import pandas as pd               # Useful to work with DataFrames and also to manage time
    from glob import glob             # To load a lot of files with same extension
    import numpy as np                # Using arrays and loading txt files
    
    # Getting path from all the files (2005-2021)
    files = glob(path+'46047*.txt')
    files.sort()
    files=np.array(files)
    years = np.arange(2015,2021,1)
    if (year_start != None) and (year_end != None): # Selecting files from the chosen years
        files = files[(years>=year_start)&(years<=year_end)]
        
    # Loading loop for the selected files
    for f in range(len(files)):
        # Loading time
        dat = np.loadtxt(files[f],skiprows=3)
    
        # Converting time vector
        timev = [str(datetime.datetime(int(dat[i,0]),int(dat[i,1]),int(dat[i,2]),int(dat[i,3]),int(dat[i,4]))) for i in np.arange(dat.shape[0])]
        timev = pd.to_datetime(timev)

        # Loading wind speed, wave height, water temperature and air temperature
        wspd = dat[:,6]
        wvht = dat[:,8]
        wtemp = dat[:,14]
        atemp = dat[:,13]
        
        # Concatenating data from different files
        if f==0:
            TIMEV,WSPD,WVHT,WTEMP,ATEMP=np.array(timev),np.array(wspd),np.array(wvht),np.array(wtemp),np.array(atemp)
        else:
            TIMEV = np.concatenate((TIMEV,np.array(timev)))
            WSPD = np.concatenate((WSPD,np.array(wspd)))
            WVHT = np.concatenate((WVHT,np.array(wvht)))
            WTEMP = np.concatenate((WTEMP,np.array(wtemp)))
            ATEMP = np.concatenate((ATEMP,np.array(atemp)))
                
    return pd.DataFrame({'wspd':WSPD,'wvht':WVHT,'wtemp':WTEMP,'atemp':ATEMP,'time':TIMEV}) # Pier object


# Writting function to perform the least squares fit to an annual cycle
def anualfit(var,t):
    from scipy.linalg import inv
    A = np.array([np.ones(t.shape[0]), np.sin(2*np.pi*t/12), np.cos(2*np.pi*t/12) ]).T
    x = np.dot(inv(np.dot(A.T, A)), np.dot(A.T,var)) # Coefficients
    fit = np.dot(A, x) # Fit
    return fit,x


# Writting function to perform the least squares fit to an annual and a semi- annual cycle
def anualsemianualfit(var,t):
    from scipy.linalg import inv
    A = np.array([np.ones(t.shape[0]),np.sin(2*np.pi*t/12), np.cos(2*np.pi*t/12),np.sin(2*np.pi*t/6), np.cos(2*np.pi*t/6) ]).T
    x = np.dot(inv(np.dot(A.T, A)), np.dot(A.T,var)) # Coefficients
    fit = np.dot(A, x) # Fit
    return fit,x

