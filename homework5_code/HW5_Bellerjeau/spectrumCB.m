function [f,a] = spectrumCB(time, data, chunk, detrend, Hanning)
%%% This function will compute a spectrum of data and return a frequency
%%% vector for plotting the spectrum
%%%     Inputs: time - time vector [n x 1]
%%%             data - data vector [n x 1]
%%%             chunk - length of chunks in indices
%%%             detrend - 1=yes, 0=no
%%%             Hanning - 1=yes, 0=no
%%%    Outputs: f - frequency vector [1 x m]
%%%             a - spectrum [1 x m]


% split into chunks
ind1 = 1;  %first index
for i=1:floor(chunk\length(data))*2  %for however many overlapping chunks will for into the data
    if ind1+chunk<length(data)  %unless we've exceeding the length of the dataset
        data1(:,i) = data(ind1:ind1+chunk);  %add a column vector to data1 with the next chunk
    end
    ind1 = ind1+chunk/2;  %step index by forward by a half-chunk
end

% frequency vector
dt = nanmean(diff(time)); %time between samples [days]
fn = 1/2/dt;  % Nyquist frequency
N= length(data1(:,1));
T = dt*N;   %length of record [days]
df = 1/T;   % fundamental frequency [cpd]
f = 0:df:fn; % frequency vector [cpd]


% compute spectrum of each chunk and average
if detrend && ~Hanning   %if the user wants detredning ubt not Hanning-ing
    for i=1:length(data1(1,:))
        data2(:,i) = data1(:,i)-mean(data1(:,i));  %just detrend
        a = fft(data2(:,i));
        amp=abs(a(1:(N+1)/2)).^2; %take half of spectrum and square
        amp = amp / N.^2; % MATLAB normalization
        amp = amp .* 2; % lost variance
        amp = amp / df; % definition of the spectrum
        A(:,i) = amp;
    end
    a = mean(A,2);  %average all spectra
end

if detrend && Hanning
    for i=1:length(data1(1,:))
        data2(:,i) = (data1(:,i)-mean(data1(:,i))).*(hann(chunk+1).*sqrt(8/3));  %detrend and Hanning
        a = fft(data2(:,i));
        amp=abs(a(1:(N+1)/2)).^2; %take half of spectrum and square
        amp = amp / N.^2; % MATLAB normalization
        amp = amp .* 2; % lost variance
        amp = amp / df; % definition of the spectrum
        A(:,i) = amp;
    end
    a = mean(A,2);  %average all spectra

end


end
