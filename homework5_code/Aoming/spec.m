function spectrum=spec(data)
%input: segmented data
%This function compute the spectrum from segmented data

dt = 10/(60*24);
fn = 1/(2*dt);
N = length(data(:,1));
T = dt*N;
df = 1/T;
f = 0:df:fn;

for i=1:length(data(1,:))
    a = fft(data(:,i));
    amp=abs(a(1:(N/2+1))).^2; %for even N
    amp = amp / N.^2;
    amp = amp .* 2;
    amp = amp / df;
    amp0(:,i) = amp;
end

spectrum = mean(amp0,2);