function [data,p] = Class10_Ben(ChunkSize)
%Class10 Plots mean spectrum
load class10_record.mat %load data
data = class10_record;
data.dnum = datevec(data.time);
NChunks = floor(length(data.time) / ChunkSize); %calc number of chunks
data.data = detrend(data.data); %detrend data
data.Chunks = reshape(data.data,ChunkSize,NChunks);
data.dt = (30/60)/24;
N=size(data.Chunks);

[data.spectra,data.f_vector,data.df] = Spectra2(data.Chunks,data.dt,'none');

[p] = Parseval(data.data,data.spectra,data.df);

% fig = figure(1);
loglog(data.f_vector,abs(data.spectra))
xlabel('Cycles per Day')
ylabel('Fourier Amplitude')
title(['Mean fft for ', num2str(NChunks), ' Chunks of ', num2str(N), ' Data Points'])
end

