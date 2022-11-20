function [err_low,err_high] = MonteCarlo_errbars(Hanning, cosine)
% This function will compute error bars on a white noise spectrum and plot
% the PDF of possible spectral levels for white noise. The user has the
% option to compute the spectra with a Hanning window or a cosine window
% Inputs: Hanning - optional Hanning windowing; 1=yes,0=no
%         cosine - optional cosine windowing; 1=yes,0=no
% Outputs: err_low - low side of error bar
%          err_high - high side of error bar

%% Generate Random Data
% generate 200 random gaussian white noise datasets
data = randn(10000,200); %each column is a dataset with 10000 points
dt = 1; %sampling rate [1 second]
N = length(data(:,1));   %length of record
time = 1:dt:N;    %time vector [seconds]

%% Compute Spectra
% compute the spectra of those 200 datasets
detrend = 1;
chunk = 500;
for i= 1:length(data(1,:))
    [f,spectra(:,i)] = spectrumCB(time,data(:,i),chunk,detrend,Hanning,cosine);
end
%each column of spectra is a spectrum, each row a frequency

if cosine && Hanning
    disp('User must select either cosine or Hanning windows, not both.')
end

%% Plot PDF
% plot PDF of all power densities obtained in 200 spectra
specvec = reshape(spectra,[50200 1]);
[NumInBin,Bins]=hist(specvec,length(f));
PDF = NumInBin./sum(NumInBin);

figure()
plot(Bins,PDF)
title('Spectrum PDF')
xlabel('\Phi_v (m/s)^2/cpd')
ylabel('Probability')


%% Calculate error bars
% for each frequency, sort the 200 power densities and find the error bars
% for 95% confidence

%preallocate
err_low = zeros(length(spectra(:,1)),1);
err_high = zeros(length(spectra(:,1)),1);

for  i=1:length(spectra(:,1))
    freqs = sort(spectra(i,:));  %sort 200 spectra at each frequency
    err_low(i) = freqs(6);  % low error bar is at 6th index
    err_high(i) = freqs(195);
end



end