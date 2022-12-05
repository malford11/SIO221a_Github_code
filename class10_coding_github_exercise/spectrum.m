%HW #7
%inputs: d = detrended data, dt = sample interval, p = chunk size
%outputs: amp = fourier amplitude, f = spectrum frequencies,
%   check_par = should be about equal to 1 if Parseval's theorem holds

function [amp, f, check_par] = spectrum(d,dt,p,window)
    N = length(d); %total number of samples
    M = N/p; %M is the number of segments
    %data = reshape(d,p,M); %reshape data into matrix of chunks
    data = buffer(d,p,p/2);
    T=dt*p; %sample period
    df=1/T; %frequency interval
    fn=1/2/dt; %nyquist frequency
    f=0:df:fn; %returned spectrum frequencies
%     window = hanning(p).*sqrt(8/3);
    data = data.*window;

    a=fft(data); %take fourier transform of data
    %normalization:
    amp=abs(a(1:p/2+1,:)).^2; % for even N
%   amp=abs(a(1:(N+1)/2).^2; % for odd N
    amp = amp / p.^2; % first correct for the MATLAB normalization
    amp = amp .* 2; %we threw out half of the spectrum; so correct for the lost variance.
    amp = amp / df; % this is then the definition of the spectrum
    data = nanmean(data,2);
    
    %check parseval's: variance of data should be equal to sum of spectrum:
    variance = nanmean(d.^2); %find variance of original detrended data
    sum_spec = nanmean(sum(amp)*df); %find sum of spectrum
    check_par = sum_spec/variance; %should be close to 1
end