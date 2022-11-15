#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def paige_spectra(data,dt,segments=1): 
    
    """
    Inputs: 
    data, 1D array of data 
    dt, sample interval [cycles per day]
    segments, integer
    
    Outputs: 
    spec = 1D array of spectra for detrended, demeaned, and overlapped data 
    freq = 1D array of frequencies 
    detrended = 2D array of detrended, demeaned, and overlapped data
    ratio = float (checking correctness of answers with Parseval's)
    """

    
    """
    Compute a 50% overlap in the data 
    """
    N = len(data) # total number of points in the data, integer
    M = segments # chunk/segments size, integer
    p = N//M # number of points in a segment, integer 
    
    data_segs = data.reshape(M,p) # reshape data into M number of segments with length p 
    
    data_ovrlp = data[p//2:-p//2].reshape(M-1,p) # reshape data into M-1 segments of length p with 50% overlap
 
    data_stack = np.vstack([data_segs,data_ovrlp]) # stack the non-overlapping and overlapping segments together
    
    
    """
    Demean and detrend the data:
    loop through every row of the overlapped and stacked data and 
    subtract the mean from every value in that row
    """
    demean = [] 

    for i in range(data_stack.shape[0]):
        data_demean = data_stack[i] - np.nanmean(data_stack[i])
        demean.append(data_demean)

    detrended = detrend(np.array(demean),axis=1)


    """
    Compute the Fourier transform of demeaned and detrended data 
    and compute the spectra 
    """    
    T = dt*p # period 
    df = 1/T # scaling factor for frequency
    fn = 1/(2*dt) # Nyquist frequency 
    
    freq = np.fft.fftfreq(p,d=dt) # frequencies array 
    freq = freq[0:int(p/2)]

    fft = np.fft.fft(detrended,axis=1) # Fourier tranform the data 
    
    amp = (abs(fft[:,0:int(p/2)])**2)/p**2
    amp = amp*2 # correct for throwing out negative values 
    specs = amp/df

    spec = np.nanmean(specs,axis=0) # take the average of the spectra in every segment
    
    """
    Verify answers with parseval's theorem  
    """
    var = detrended.var() # variance of the data 
    pars = np.sum(spec)*df # integral of the spectra multiplied by the fundamental frequency 
    ratio = var/pars # ratio of the variance and integral of spectra 
    
    return spec,freq,detrended,ratio

