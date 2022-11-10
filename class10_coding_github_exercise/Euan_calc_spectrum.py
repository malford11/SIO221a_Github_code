def Euan_calc_spectrum(data, sample_interval, chunk_size=None, no_chunks=None, data_units=None, time_units=None):
    
    import numpy as np
    from scipy import signal
    
    """
    CALCULATES SPECTRUM OF A TIMESERIES
    
    Inputs
    ------
    data: array
        1D array of data
    sample_interval: float
        Sample interval in units of choice
    chunk_size = int
        Length of each segment
    no_chunks = int
        Number of chunks. If both no_chunks and chunk_size are specified then no_chunks takes priority
        
    Outputs
    -------
    f: array
        Array of frequency values
    absX: array
        Array of power spectrum values
    
    """
    
    length = len(data)
            
    if no_chunks:
        chunk_size = int(np.floor(length/no_chunks))
    elif chunk_size:
        no_chunks = int(np.floor(length/chunk_size))
    
    # Print the data length, number of chunks and size of each chunk
    print('The length of the data is', length)
    print('The number of segments in the data is', no_chunks)
    print('The size of each segments is', chunk_size)

    ending = int(chunk_size*no_chunks)
    data_snipped = data[:ending]
    
    # Divide the data into chunks
    data_chunked = np.zeros((chunk_size, no_chunks))
    for i in range(no_chunks):
        data_chunked[:, i:i+1] = np.expand_dims(data[int((i)*chunk_size):int((i+1)*chunk_size)], axis=1)
    
    # Detrend the data to ensure Parseval's theorm holds
    data_detrended = signal.detrend(data_chunked, axis=0)
    
    # N is the number of points in each chunk
    # M is the number of chunks
    N = np.shape(data_chunked)[0]
    M = np.shape(data_chunked)[1]
    
    dt = sample_interval
    T = dt*N      
    df = 1/T         
    fn = 1/(2*dt)

    # Find the frequency vector ranging from 0 to the Niquist frequency
    f = np.arange(0, fn, df)
    
    # Fourier transform the timeseries data
    X = np.fft.fft(data_detrended, axis=0)
        
    absX = abs(X[0:int((N+1)/2),:])**2

    absX = absX*2 # We threw out half of the spectrum; so correct for the lost variance
    absX = absX/N**2 # First correct for the normalization
    absX = absX/df # This is then the definition of the spectrum
    
    # Test Parseval's theorem
    variance = np.mean(np.mean(data_detrended**2, axis=1))
    sum_spec = np.sum(np.mean(absX, axis=1))*df
    
    print('\nThe ratio between the sum of the spectrum and the variance is', np.round(sum_spec/variance,5))
    
    if abs(sum_spec/variance-1) < 0.005:
        print("Parseval's theorem holds\n")
    else:
        print("Parseval's theorem does not hold\n")
      
    return f, absX