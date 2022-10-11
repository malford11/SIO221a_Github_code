def GetPierDataMY(path,year_start=None,year_end=None):
    # This function returns a multi-year temperature, pressure and time pier records 
    # (http://sccoos.org/thredds/catalog/autoss/catalog.html)
    
    # The path must be a directory containing all the pier data from 2005-2021
    # year_start is the start year of interest (must be between 2005 and 2021)
    # year_end is the end year of interest (must be between 2005 and 2021)
    # Default is year_start and year_start equal to None (all the 2005-2021 record will be loaded)
    # CDL SIO211A 10/08/2022

    # Important packages
    import netCDF4 as nc              # Read nc files
    import datetime                   # Manage time vectors
    import pandas as pd               # Useful to work with DataFrames and also to manage time
    from glob import glob             # To load a lot of files with same extension
    import numpy as np                # Using arrays
    
    # Getting path from all the files (2005-2021)
    files = glob(path+'scripps*.nc')
    files.sort()
    files=np.array(files)
    years = np.arange(2005,2022,1)
    if (year_start != None) and (year_end != None): # Selecting files from the chosen years
        files = files[(years>=year_start)&(years<=year_end)]
        
    # Loading loop for the selected files
    for f in range(len(files)):
        # Loading time
        dat = nc.Dataset(files[f])
        time = dat['time'][:]
    
        # Converting time vector
        timev = [str(datetime.datetime(1970,1,1)+datetime.timedelta(seconds=int(d))) for d in time]
        timev = pd.to_datetime(timev)
        # Loading temp and pressure
        temp = dat['temperature']
        pres = dat['pressure']
        
        # Concatenating data from different files
        if f==0:
                   TIMEV,TEMP,PRES=np.array(timev),np.array(temp),np.array(pres)
        else:
                   TIMEV = np.concatenate((TIMEV,np.array(timev)))
                   TEMP = np.concatenate((TEMP,np.array(temp)))
                   PRES = np.concatenate((PRES,np.array(pres)))
                
    return pd.DataFrame({'temp':TEMP,'pres':PRES,'time':TIMEV}) # Pier object



# Creating function to compute the mean and standard error for temperature and pressure records from a particular year
def mean_stderror(temp,pres,timev,year):
    # Returns the mean and standard error for temperature and pressure records from a particular year
    # CDL 10/09/2022
    
    import numpy as np
    import datetime
    
    temp = temp[timev.year==year] # Slicing vector for the selected year
    pres = pres[timev.year==year] # Slicing vector for the selected year
    
    # Computing means
    tmean = np.nanmean(temp)
    pmean = np.nanmean(pres)

    # Computing standard deviations
    tstd = np.nanstd(temp)
    pstd = np.nanstd(pres)

    # N
    Nt = temp[~np.isnan(temp)].shape[0] # Size of temp matrix without NaNs 
    Np = pres[~np.isnan(pres)].shape[0] # Size of temp matrix without NaNs 

    # Computing standard error of the mean
    tstderr = tstd/(Nt**(1/2))
    pstderr = pstd/(Np**(1/2))

    return tmean,pmean,tstderr,pstderr


# Creating function to compute the variance and its standard error for temperature and pressure records from a particular year
def var_stderror(temp,pres,timev,year):
    # Returns the variance and its standard error for temperature and pressure records from a particular year
    # CDL 10/09/2022
    
    import numpy as np
    import datetime
    
    temp = temp[timev.year==year] # Slicing vector for the selected year
    pres = pres[timev.year==year] # Slicing vector for the selected year
    
    # Computing variances
    tvar = np.nanstd(temp)**2
    pvar = np.nanstd(pres)**2

    # N
    Nt = temp[~np.isnan(temp)].shape[0] # Size of temp matrix without NaNs 
    Np = pres[~np.isnan(pres)].shape[0] # Size of temp matrix without NaNs 

    # Computing standard error of the mean
    tstderr = tvar*(2/(Nt-1))**(1/2)
    pstderr = pvar*(2/(Np-1))**(1/2)

    return tvar,pvar,tstderr,pstderr


# Using function from lecture to compute histogram

def compute_histogram(variable, bin_max, bin_min, dbin, pdf=False):
    
    """ Computes 1D histogram or probability density for a given variable.
        
    Keyword arguments:
    variable -- 1D array.
    bin_max -- maximum value for bins
    bin_min -- minimum value for bins
    dbin -- bin size
    pdf -- (default False)
    
    Returns:
    bins -- histogram bins
    counts -- either counts or probability density
        
    """
    bins = np.arange(bin_min, bin_max, dbin)
    count = []
    for i in range(len(bins)):
        ind = (variable>bins[i] - dbin/2) & (variable<=bins[i]+dbin/2)
        count.append(ind.sum())
    count = np.array(count)
    if pdf:
        norm_hist = count/count.sum()/dbin
        assert np.allclose(norm_hist.sum()*dbin, 1.0), "PDF doesn't sum to 1"
    
        return bins, norm_hist
    else:
        return bins, count   

# Creating a function to compute the likelihood of a valuer being higher than the mean + level*sigma
def likelihood(var,varmax,varmin,dbin,level,lower=False,upper=True):
    
    bins,pdf = compute_histogram(var,varmax,varmin,dbin,pdf=True)
    cdf = np.cumsum(pdf)
    varmean = np.nanmean(var)
    varstd = np.nanstd(var)
    cut = varmean+level*varstd
    return 1 - cdf[np.where(bins>=cut)[0][0]]


# Defining function to compute gaussian function
def gaus(var,bins):
    sig = np.nanstd(var)
    mu = np.nanmean(var)
    return (1.0/sig/np.sqrt(2*np.pi)) * (np.exp(-(bins-mu)**2 / (2*sig**2)))


