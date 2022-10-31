
def seg(var,nchunk):
	# 1. ensure nchunk divides evenly into length of var
	assert var.shape[0]%nchunk == 0
	# 2. reshape so that we have segsize * length/segzize
	data = np.reshape(var, nchunk,len(var)/nchunk)
	# 3. return segmented data
	return data
from scipy.signal import detrend


def fftnorm(fft,N,tstep):
	ffts = abs(fft)/N**2 
	ffts = ffts**2/(1/tstep) 
	return ffts

def fftseg(var,nchunk,tstep):
	segs = seg(var,nchunk)
	freq = np.fft.fftshift(np.fft.fftfreq(segs.shape[1],d=tstep))
	fft = np.fft.fftshift(np.fft.fft(segs),axes=1)

	idx = np.where(freq>0)[0]
	freq = freq[idx]
	fft = fft[:,idx]
	fft = fftnorm(fft,segs.shape[1],tstep)
	
	return freq,fft
