function [f,a,Parseval] = spectrumCB_demo(time, data, chunk)

% split into chunks
ind1 = 1;  %first index
for i=1:floor(chunk\length(data))*2  %for however many overlapping chunks will for into the data
    if ind1+chunk<length(data)  %unless we've exceeding the length of the dataset
        data1(:,i) = data(ind1:ind1+chunk);  %add a column vector to data1 with the next chunk
    end
    ind1 = ind1+chunk/2;  %step index by forward by a half-chunk
end

%frequencies
dt = nanmean(diff(time)); %time between samples [days]
fn = 1/2/dt;  % Nyquist frequency
N= length(data1(:,1));
T = dt*N;   %length of record [days]
df = 1/T;   % fundamental frequency [cpd]
f = 0:df:fn; % frequency vector [cpd]


%compute spectrum of each chunk and average
for i=1:length(data1(1,:))
    data2(:,i) = detrend(data1(:,i));
    a = fft(data2(:,i));
    amp=abs(a(1:(N+1)/2)).^2; %take half of spectrum and square
    amp = amp / N.^2; % MATLAB normalization
    amp = amp .* 2; % lost variance
    amp = amp / df; % definition of the spectrum
    A(:,i) = amp;
end
a = mean(A,2);

%check Parseval's theorem
variance=std(data)^2;
int_spec = trapz(f,a);
%int_spec = sum(a)*df;
Parseval=int_spec/variance;

end
