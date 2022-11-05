def compspectrum(data,tday,tstep,demean=True,detren=False,hanning=False):
    # data is an array with the data
    # tday is the time vector in days
    # tstep should be in days

    import numpy as np
    import xarray as xr
    from scipy.signal import detrend

    # Segmenting
    segs=np.arange(0,int((tday[-1]/60 -1)//0.5),1)
    tdays = np.arange(0,60+tstepd,tstepd)
    blank = np.ones((segs.shape[0],tdays.shape[0]))*np.nan

    # Creating xarray with dimensions segment and time 
    xr1 = xr.Dataset({"variable": (("Segment", "Time"),blank.copy())},coords={"Segment":segs,"Time": tdays})
    n30=np.where(tdays==30)[0][0]
    for i in np.array(xr1.Segment):
        xr1.variable[i] = np.array(data[i*n30:((i+2)*n30)+1])

    data = xr1.variable.copy()
    
    # Removing the mean
    if demean == True:
        data = data-data.mean(dim='Time')
        
    # Computing frequencies
    freq = np.fft.fftshift(np.fft.fftfreq(data.shape[1],d=tstep))
    # Indexes of positive freqs
    idx = np.where(freq>0)[0]
    freq=freq[idx]
    
    if detren==False and hanning==False:
        ffts = np.fft.fftshift(np.fft.fft(data),axes=1)[:,idx]
    if detren==True and hanning ==False:
        ffts = np.fft.fftshift(np.fft.fft(detrend(data)),axes=1)[:,idx]
    if detren==True and hanning == True:
        hann = np.hanning(data.shape[1])
        fac = np.sqrt(np.nanmean(hann**2)) 
        hann = hann/fac # Normalizing to get <hann**2>=1
        Hann,_ = np.meshgrid(hann,data.Segment)
        ffts = np.fft.fftshift(np.fft.fft(detrend(data)*Hann),axes=1)[:,idx]

    ffts = abs(ffts)/data.shape[1]**2 
    ffts = ffts**2/(1/tstep) 

    return np.nanmean(ffts,axis=0),freq,segs.shape[0]


def seg60day(df,tday,tstepd):
    segs=np.arange(0,int((tday[-1]/60 -1)//0.5),1)
    tdays = np.arange(0,60+tstepd,tstepd)
    blank = np.ones((segs.shape[0],tdays.shape[0]))*np.nan
    import xarray as xr
    xr1 = xr.Dataset({"wspd": (("Segment", "Time"),blank.copy()),"uwnd": (("Segment", "Time"),blank.copy()),"vwnd": (("Segment", "Time"),blank.copy())},coords={"Segment":segs,"Time": tdays})
    n30=np.where(tdays==30)[0][0]
    for i in np.array(xr1.Segment):
        xr1.wspd[i] = np.array(df.wspd[i*n30:((i+2)*n30)+1])
        xr1.uwnd[i] = np.array(df.uwnd[i*n30:((i+2)*n30)+1])
        xr1.vwnd[i] = np.array(df.vwnd[i*n30:((i+2)*n30)+1])
        
    return xr1

def subplot_time(timev,var,varstr,ax,c,labels=None):
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    if labels!=None:
        plt.plot(timev,var,label=labels,color=c)
        plt.legend(loc='best')
    else:
        plt.plot(timev,var,color=c)
    plt.xlim(timev[0],timev[-1])
    plt.xlabel('Time [MM-DD-YY]',fontsize=12,weight='bold')
    plt.ylabel(varstr,fontsize=12,weight='bold')
    # Setting datetick format
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d-%Y'))
    for label in ax.get_xticklabels(which='major'):
        label.set(rotation=20, horizontalalignment='right') 


def GetTAOPData(file):

    # Important packages
    import netCDF4 as nc              # Read nc files
    import datetime                   # Manage time vectors
    import pandas as pd               # Useful to work with DataFrames and also to manage time
    import numpy as np                # Using arrays
    
    # Loading time
    dat = nc.Dataset(file)
    time = dat['TIME'][:]
    # Converting time vector
    timev = [str(datetime.datetime(1950,1,1)+datetime.timedelta(days=d)) for d in time]
    timev = pd.to_datetime(timev)
    # Loading wind data
    wspd = dat['WSPD'][:,0]
    uwnd = dat['UWND'][:,0]
    vwnd = dat['VWND'][:,0]
        
    df = pd.DataFrame({'wspd':wspd.data,'uwnd':uwnd.data,'vwnd':vwnd.data})
    df.index = timev # Indexing DataFrame with timev
                
    return df 


