function [x,y] = MidtermFunc(data,Hanning,Detrend,Chunked)
%data inputs:
%NonChunked is 1xn time series.
%Chunked is 86400xn chunked segments
%Options: WSPD, VWND, UWND
%years: Options1: 2015, Options2: 2016, Options3: 2017, Options: Full
%Time Series.
%Example: Chunks.WSPD1: 2015 Chunked Windspeed
%Example2: NonChunked.VWND: Full 1xn nonchunked time series
%0 or 1 to apply hanning and detrend options.
%if using chunked data must put 1 for 'Chunked' option.
%x outputs a frq vector
%y outputs the corresponding spectra

x = frqvector(data,10);
y = Spectra(data,Detrend,Hanning,Chunked);

end

