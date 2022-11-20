#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import xarray as xr
from scipy import signal
from numpy.fft import fft, ifft, fftfreq

def spectra(da,sr,method,npts,overlap=True):
    
    """
    Function to calculate power spectral density for a given data array and sampling rate.
    
    --- INPUTS ---
    da: input array of data, must be in list format (e.g. [data], or if multiple [data1,data2,...])
    sr: sampling rate of data  (units are implied here)
    method: raw, detrended, or hanning window
    npts: number of points per segments
    overlap: True if 50% overlapping segments are desired, False if not
    
    --- OUTPUTS ---
    
    freq: array with frequencies over which spectrum was calculated
    pwr: array of spectra for each segment - need to take the mean over axis 0 to get mean spectrum
    da_segments: array of input data after it has been segmented (can see how many segments were created from this)
    """
    
    seg_arr_list = [] # create empty list for segmented arrays
    nsegments_total = 0

    for i in da: # loop through list of data arrays 
        
        if overlap is False:
            
            nsegments = int(len(i)/(npts)) 

            seg_arr = np.zeros((nsegments,npts)) # empty array of shape (nsegments,npts) 

            for k,j in enumerate(range(0,nsegments)): 
                seg = i[int(0+j*(npts)):int(npts+j*(npts))]
                seg_arr[k,:] = seg 

        else:
            
            nsegments = int(len(i)/(npts/2)) - 1 # need to divide npts by 2 to account for 50% overlap

            seg_arr = np.zeros((nsegments,npts)) # empty array of shape (nsegments,npts) 

            for k,j in enumerate(range(0,nsegments)): 
                seg = i[int(0+j*(npts/2)):int(npts+j*(npts/2))] # 50% overlap formula
                seg_arr[k,:] = seg 
            
        seg_arr_list.append(seg_arr) # append each segment to array of segments
        nsegments_total = nsegments_total+nsegments
    
    da_segments = np.vstack(seg_arr_list) # combine all segment arrays into one big segment array (in case multiple years data used)
    
    N = npts # Number of data points
    T = N/sr # Period of segment

    fn = sr/2 # Nyquist frequency
    fo = 1/T # Fundamental frequency
    df = fo # Frequency bin width 

    freq = np.arange(0,fn,df) # frequency space
        
    if method == 'raw':
        
        da_processed = da_segments
    
    elif method == 'detrend':
        
        da_processed = signal.detrend(da_segments) # detrend and demean
        
    elif method == 'hanning':
        
        da_processed = signal.detrend(da_segments)
        w = np.hanning(N) # hann window
        w_norm = (w**2/np.mean(w**2))**0.5 # normalizing hann window so that its squared mean is 1
        da_processed = w_norm*da_processed # data multiplied by normalized hann window
     
    elif method == 'cosine':

        da_processed = signal.detrend(da_segments)
        time = np.arange(-T/2,T/2,(1/sr))
        c = np.cos(np.pi*time/T)
        c_norm = (c**2/np.mean(c**2))**0.5
        da_processed = c_norm*da_processed # data multiplied by normalized hann window
        
    X = fft(da_processed) # fourier transform of each data segment

    ampsq = abs(X[:,0:int((N+1)/2)])**2 # finding spectral energy, only need one side of spectrum
    
    ampsq = ampsq/(N**2) # normalizing FFT

    ampsq = ampsq*2 # correct for variance lost by using half of spectrum
    
    pwr = ampsq/df # spectral power is spectral energy divided by frequency bin width
 
    ### Parseval's theorem check ###
    time_array = signal.detrend(da_segments)
    
    time_pars =  np.mean(np.mean(time_array**2,axis=0)) # time domain
    freq_pars = np.sum(np.mean(pwr,axis=0))*df # frequency domain
    ratio = time_pars/freq_pars
    
    if abs(ratio-1)<0.01:
        print ('The ratio calculated is {:.3f}.'.format(ratio))
        print ("Parseval's theorem is verified.")
        
    else:
        print ('The ratio calculated is {:.3f}.'.format(ratio))
        print ("Parseval's theorem does not hold.")
        
    return freq,pwr,da_segments