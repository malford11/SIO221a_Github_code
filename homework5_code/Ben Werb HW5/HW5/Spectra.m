function [y] = Spectra(data,Detrend,Hanning,Chunked)
%Output the Spectra to Plot with Hanning Window Applied
n = length(data);
x = size(data);
x = x(2);
y = data;
if Detrend>0
    y = detrend(y,1);
end
if Hanning>0 && Chunked>0
    hanningwindow = repmat(hann(n),1,x) * (8/3);
    y = y .* hanningwindow;
elseif Hanning>0 && Chunked<0
    hanningwindow = hann(n) * (8/3);
    y = y .* hanningwindow;
end
y = fftshift(fft(y));
if Chunked>0
    y = nanmean(y,2);
end
y = abs(y);
end