function [f, spectrum, ratio]=spec_lec10(data,time, seg_size)
%input: data, time, segment size
%This function compute the spectrum of data with a specified segment size

data = detrend(data); %remove mean
seg_num = floor(length(data)/seg_size); %number of segments
data_seg = reshape(data, seg_size, seg_num); %reshape

dt = mean(diff(time)); %time interval
fn = 1/(2*dt); %Nyquist frequency
N = length(data_seg(:,1)); %segment length
T = dt*N; %time span
df = 1/T; %frequency resolution
f = 0:df:fn; %frequency vector

% Calculate spectrum for each segment
for i=1:length(data_seg(1,:))
    a = fft(data_seg(:,i));
    amp=abs(a(1:(N/2+1))).^2; %for even N
    amp = amp / N.^2;
    amp = amp .* 2;
    amp = amp / df;
    amp0(:,i) = amp;
end

% Mean of all segments spectrum
spectrum = mean(amp0,2);

% Parseval's Theorem
variance = mean(data_seg.^2);
sum_spec = sum(amp0)*df;
ratio = sum_spec / variance;
end
