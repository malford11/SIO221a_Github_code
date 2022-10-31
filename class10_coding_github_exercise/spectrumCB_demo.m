function [f,a] = spectrum(time, data, chunk)

chunk = 1600; %hardcoding for this class demo ;)

%split into chunks
ind1 = 1;
for i=1:10
    data1(:,i) = data(ind1:ind1+chunk);
    ind1 = ind1+chunk/2;
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
    data2(:,i) = data1(:,i)-mean(data1(:,i));
    a = fft(data2(:,i));
    amp=abs(a(1:(N+1)/2)).^2; %take half of spectrum and square
    amp = amp / N.^2; % MATLAB normalization
    amp = amp .* 2; % lost variance
    amp = amp / df; % definition of the spectrum
    A(:,i) = amp;
end
a = mean(A,2);


end
