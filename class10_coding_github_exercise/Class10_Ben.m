function [data,p] = Class10_Ben(ChunkSize)
%Class10 Plots mean spectrum
load class10_record.mat %load data
data = class10_record; %make a nice struct
data.dnum = datevec(data.time); %easy to read dates
NChunks = floor(length(data.time) / ChunkSize); %calc number of chunks
data.data = detrend(data.data); %detrend data
data.Chunks = reshape(data.data,ChunkSize,NChunks); %chunkify
data.dt = (30/60)/24; %convert dt=30 min into CPD
N=size(data.Chunks);

[data.spectra,data.f_vector,data.df] = Spectra2(data.Chunks,data.dt,'none');
%calculate spectra
[p] = Parseval(data.data,data.spectra,data.df);
%check Parseval's theorem
fig = figure(1);
loglog(data.f_vector,abs(data.spectra))
xlabel('Cycles per Day')
ylabel('Fourier Amplitude')
title(['Mean fft for ', num2str(NChunks), ' Chunks of ', num2str(N), ' Data Points'])
end

