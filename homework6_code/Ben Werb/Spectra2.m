function [spectra,frqvector,df,fn] = Spectra2(data,dt,Window)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
S = size(data);
if S(2)>1
    data=detrend(data,1);
elseif S(2)==1
    data=detrend(data);
end
N=S(1);

if ismember(Window,'None')
    data=data;
elseif ismember(Window,'Hanning')
    hanningwindow = repmat(hann(N),1,S(2)) * (8/3);
    data = data .* hanningwindow;
elseif ismember(Window,'Cosine')
    %(cos(πt/T ), where −T /2 < t < T /2).
    t = linspace(-0.5,0.5,500);
    T = 1; %period is 1
    cos_window = cos(pi*t/T); %I used Okun's HW on Github as a guide for
    %this. I could not figure out how to get the Cos window working
    %correctly.
    cos_factor = sqrt(sum(cos_window.^2)/500);
    cos_window = cos_window'*(1/cos_factor);
    CosWindow = repmat(cos_window,1,S(2));
    data = data .* CosWindow;
%     CosWindow = repmat(tukeywin(N,.25),1,S(2)) * (8/3); 
%     data = data .* CosWindow;
end

if S(2)==1
    a=fft(data);
    amp=abs(a(1:N/2+1)).^2; % for even N
    %  amp=abs(a(1:(N+1)/2).^2; % for odd N
    dt=(dt/60)/24;
    T=dt*N;
    df=1/T;
    fn=(1/2)/dt;
    frqvector=0:df:fn; %frequency vector, cpd, goes from 0 to Nyquist.
    amp = amp / N.^2; % first correct for the MATLAB normalization
    amp = amp .* 2; %we threw out half of the spectrum; so correct for the lost variance.
    spectra = amp / df; % this is then the definition of the spectrum
elseif S(2)>1
    p = S(1);
    a=fft(data);  % this computes the fft for each column
    amp=abs(a(1:p/2+1,:)).^2;
    dt=(dt/60)/24;
    T=dt*p;
    df=1/T;
    fn=1/2/dt;
    frqvector=0:df:fn; %frequency vector, cpd, goes from 0 to Nyquist.
    % Normalize as above
    amp = amp / p.^2; % first correct for the MATLAB normalization
    amp = amp .* 2; %we threw out half of the spectrum; so correct for the lost variance.
    amp = amp / df; % this is then the definition of the spectrum
    spectra = nanmean(amp,2);
end
end

