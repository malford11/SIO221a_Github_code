def compspectrum_caique(var,tstep,nchunk=1,detren=False):
    # data is an array with the data
    # tstep should be in days
    # nchunk is the number of chunks
    # detrend needs to be True to demean and detrend the data
    # Caique Luko

    import numpy as np
    from scipy.signal import detrend

    # 1. ensure nchunk divides evenly into length of var
    assert var.shape[0]%nchunk == 0

    # 2. reshape so that we have seg x length/seg
    data = var.copy().reshape(int(nchunk),int(len(var)/nchunk))

    # 3. Removing the mean and detrending
    if detren == True:
        data = detrend(data,axis=1)

    # 4. Computing frequencies
    freq = np.fft.fftfreq(data.shape[1],d=tstep)

    # 5. Computing Fourier Transform
    ffts = np.fft.fft(data,axis=1)

    # 6. Taking only positive frequencies
    N = data.shape[1]
    freq = freq[1:int(N/2)]
    ffts = ffts[:,1:int(N/2)]
    
    # 7. Spectrum
    T = N*tstep           # Total period
    df = 1./T             # Fundamental frequency
    ffts = (abs(ffts)**2) # Spectral energy
    ffts = ffts/(N**2)    # Normalizing
    ffts = 2*ffts         # Correcting for taking only half of the Fourier coefficients
    spec = ffts/df        # Spectral power 
    
    
    # 8. Parseval
    vari = np.mean(np.std(detrend(data,axis=1),axis=1)**2)   # Variance of detrended original series
    pars = np.mean(np.sum(spec,axis=1)*df)                   # Integral of the spectrum in frequency
    ratio = vari/pars
    
    spec = np.mean(spec,axis=0) # Mean of all segments

    return spec,freq,data.shape[0],ratio # spectrum, frequency, number of segments, parseval ratio
