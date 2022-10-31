function Class10(ChunkSize)
%Class10 Plots mean spectrum
load class10_record.mat %load data
df = class10_record;
df.dnum = datevec(df.time);
NChunks = floor(length(df.time) / ChunkSize); %calc number of chunks
df.data = detrend(df.data); %detrend data
df.Chunks = reshape(df.data,ChunkSize,NChunks);
N=size(df.Chunks);
N=N(1); %this is the same as chunk size
dt=(30/60)/24; %in cycles per day
T=dt*N;
dfft=1/T;
fn=1/2/dt;
f_vector=-fn:dfft:(fn-dfft);
df.fft = fftshift(fft(df.Chunks));
df.fft_norm = df.fft / (NChunks^2) * ChunkSize; %normalize
df.meanfft_norm = mean(df.fft_norm,2); %mean

% fig = figure(1);
plot(f_vector,abs(df.meanfft_norm))
xlabel('Cycles per Day')
ylabel('Fourier Amplitude')
title(['Mean fft for ', num2str(NChunks), ' Chunks of ', num2str(N), ' Data Points'])
end

