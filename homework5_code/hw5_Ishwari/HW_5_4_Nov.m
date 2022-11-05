%SIOC 221A Data Analysis
%Homework 5
%Ishwari %November 4, 2022

% I certify that this represents my own work and that I have not worked
% with classmates or other individuals to complete this assignment

%% Preliminary assessment of the data

%local path %hardcoded
filePath = 'C:\Users\Ishwari\OneDrive\Desktop\SIOC_221A\Data\';

%2015
fileName2015 = 'OS_T8S110W_DM134A-20150425_D_WIND_10min.nc';
fileIn2015 = [filePath fileName2015];

%2016
fileName2016 = 'OS_T8S110W_DM183A-20160321_D_WIND_10min.nc';
fileIn2016 = [filePath fileName2016];

%2017
fileName2017 = 'OS_T8S110W_DM231A-20170606_D_WIND_10min.nc';
fileIn2017 = [filePath fileName2017];


%%% reading in data %%%
wspd2015 = ncread(fileIn2015, 'WSPD');
uwnd2015 = ncread(fileIn2015, 'UWND');
vwnd2015 = ncread(fileIn2015, 'VWND');
time2015 = ncread(fileIn2015, 'TIME');

wspd2016 = ncread(fileIn2016, 'WSPD');
uwnd2016 = ncread(fileIn2016, 'UWND');
vwnd2016 = ncread(fileIn2016, 'VWND');
time2016 = ncread(fileIn2016, 'TIME');

wspd2017 = ncread(fileIn2017, 'WSPD');
uwnd2017 = ncread(fileIn2017, 'UWND');
vwnd2017 = ncread(fileIn2017, 'VWND');
time2017 = ncread(fileIn2017, 'TIME');

%converting to datetime
time2015 = datetime(time2015 + datenum('1950-01-01 00:00:00'), 'ConvertFrom', 'datenum');
time2016 = datetime(time2016 + datenum('1950-01-01 00:00:00'), 'ConvertFrom', 'datenum');
time2017 = datetime(time2017 + datenum('1950-01-01 00:00:00'), 'ConvertFrom', 'datenum');

%spacing between time data
%plot(diff(time2015))
%plot(diff(vertcat(time2015,time2016)))
%% plotting raw data
fig=figure;
subplot(3,1,1)
plot(time2015, wspd2015)
legend('wind speed')
subplot(3,1,2)
plot(time2015, uwnd2015)
legend('u wind')
subplot(3,1,3)
plot(time2015, vwnd2015);
legend('v wind')
xlabel('time')
h=axes(fig,'visible','off'); 
h.YLabel.Visible = 'on';
h.Title.Visible = 'on';
ylabel(h, 'velocity in m/s')
title('Raw data plots for 2015')

%% removing gaps and interpolating

%InterpolateData replaces NaN values with value from MATLAB interp1
%2015
wspd2015 = InterpolateData(wspd2015, time2015); 
uwnd2015 = InterpolateData(uwnd2015, time2015);
vwnd2015 = InterpolateData(vwnd2015, time2015);

%2016
wspd2016 = InterpolateData(wspd2016, time2016);
uwnd2016 = InterpolateData(uwnd2016, time2016);
vwnd2016 = InterpolateData(vwnd2016, time2016);

%2017
wspd2017 = InterpolateData(wspd2017, time2017);
uwnd2017 = InterpolateData(uwnd2017, time2017);
vwnd2017 = InterpolateData(vwnd2017, time2017);


%% Segmenting data

L = 60;  %length of chunk in days
a = 0.5; %overlap fraction

%function Nchunks returns an (n_points_per_chunk *n_chunks) matrix containing chunks in columns
%function chunkify returns a vector [n_chunks, n_points_per_chunk] for given overlap and segment length

%%% this is basically what chunkify.m does: %%%
% nchunks = floor(N(1));
% npoints = N(2);
% del = npoints*(1-a);
% A = [];
% for j = 1:nchunks
%    if  floor((j-1)*del)+ npoints > length(wspd2015)
%        A(1: npoints, j) = horzcat(wspd2015(1 + floor((j-1)*del): length(wspd2015)), zeros(1,floor((j-1)*del)+ npoints - length(wspd2015)))';
%    else
%        A(1: npoints, j) = wspd2015(1 + floor((j-1)*del): floor((j-1)*del)+ npoints)';
%    end    
% end 


%2015
N = Nchunks(a, L, time2015); 
Cwspd2015 = chunkify(a, N, wspd2015);
Cuwnd2015 = chunkify(a, N, uwnd2015);
Cvwnd2015 = chunkify(a, N, vwnd2015);
%%% chunkify will not work on datetime so convert to datenum %%%
Ctime2015 = chunkify(a, N, datenum(time2015) - min(datenum(time2015))); 

%2016
N = Nchunks(a, L, time2016);
Cwspd2016 = chunkify(a, N, wspd2016);
Cuwnd2016 = chunkify(a, N, uwnd2016);
Cvwnd2016 = chunkify(a, N, vwnd2016);
Ctime2016 = chunkify(a, N, datenum(time2016) - min(datenum(time2016)));

%2017
N = Nchunks(a, L, time2017);
Cwspd2017 = chunkify(a, N, wspd2017);
Cuwnd2017 = chunkify(a, N, uwnd2017);
Cvwnd2017 = chunkify(a, N, vwnd2017);
Ctime2017 = chunkify(a, N, datenum(time2017) - min(datenum(time2017)));

df = 1/L;      %fundamental frequency in cpd
fn = 1*6*24/2; %Nyquist freq in cpd
f = 0:df:fn;   %frequency vector

%% The 3 spectra for 2015 wind speed

%%% without detrending or windowing %%%
X = fft(Cwspd2015);
amp_X = abs(X).^2;
amp_X = (2/(size(X,1)^2))*amp_X;
amp_X = amp_X(1:floor(size(X,1)/2)+1,:);
amp_X = amp_X/df;


%%% with detrend %%%
Y = fft(detrend(Cwspd2015));
amp_Y = abs(Y).^2;                       %calculating amplitude 
amp_Y = (2/(size(Y,1)^2))*amp_Y;         %correcting for MATLAB normalisation
amp_Y = amp_Y(1:floor(size(Y,1)/2)+1,:); %taking only the first half of the frequency values
amp_Y = amp_Y/df;                        %spectrum



%detrend and Hanning window
hanwin = cos(pi*Ctime2015/ L).^2;         %defining the Hanning window function
A = Cwspd2015;
A(1:end,:) = hanwin(1:end,:).*A(1:end,:); %convolving the window with data
Z = fft(detrend(A));
amp_Z = abs(Z).^2;
amp_Z = (2/(size(Z,1)^2))*amp_Z;
amp_Z = amp_Z(1:floor(size(Z,1)/2)+1,:);
amp_Z = amp_Z/df;

%%% plotting the 3 types together %%%
loglog(f, nanmean(amp_X,2),'r', f,nanmean(amp_Y,2),'b',f,nanmean(amp_Z,2),'k');
legend('raw data', 'detrend', 'window')
title('Spectrum for 2015 wind speed')
xlabel('frequency in cpd')
ylabel('spectral power in m^2/s^2 per cpd')

%% Uncertainty estimates


nu=2*size(Z,2);                   %degrees of freedom
err_low = nu/chi2inv(1-.05/2,nu); %calculating range using chi2inv function
err_high = nu/chi2inv(.05/2,nu);

%%% plotting E(f) with uncertainty for some median f value %%%

%%% windowed data %%%
fig = figure;
subplot(3,1,1)
exp = 2*sum(abs(Z(1:floor(size(Z,1)/2)+1,:)).^2,2)/size(Z,2);  %summing squares over f, by definition of E(f)
exp = exp(2:end);        %removing the expectation for zero frequency as it's too high
semilogy(f(2:end),exp)   %plotting E(f)
hold on 
%calculating and plotting uncertainty at N/2
semilogy([f(floor(length(f)/2)) f(floor(length(f)/2))], [err_low err_high]*exp(floor(length(f)/2)), 'LineWidth', 3)
xlim([0 70])
legend('detrended+windowed', 'uncertainty')

%%% detrended data %%%
subplot(3,1,2)
exp = 2*sum(abs(Y(1:floor(size(Y,1)/2)+1,:)).^2,2)/size(Y,2);
exp = exp(2:end);
semilogy(f(2:end),exp)
hold on 
semilogy([f(floor(length(f)/2)) f(floor(length(f)/2))], [err_low err_high]*exp(floor(length(f)/2)), 'LineWidth', 3)
xlim([0 70])
legend('detrended', 'uncertainty')

%%% raw data %%%
subplot(3,1,3)
exp = 2*sum(abs(X(1:floor(size(X,1)/2)+1,:)).^2,2)/size(X,2);
exp = exp(2:end);
semilogy(f(2:end),exp)
hold on 
semilogy([f(floor(length(f)/2)) f(floor(length(f)/2))], [err_low err_high]*exp(floor(length(f)/2)), 'LineWidth', 3)
xlim([0 70])
xlabel('frequency in cpd')
legend('raw data', 'uncertainty')

h=axes(fig,'visible','off'); 
h.YLabel.Visible = 'on';
h.Title.Visible = 'on';
ylabel(h, 'Expectation value of chi squared variable')
title('Uncertainty estimates at 95% confidence level')


%% Using multiyear data

Cwspd = horzcat(Cwspd2015, Cwspd2016, Cwspd2017); %creating single wspd array
Ctime = horzcat(Ctime2015, Ctime2016, Ctime2017); %creating single time array

hanwin = cos(pi*Ctime/ L).^2;             %defining Hanning window
A = Cwspd; 
A(1:end,:) = hanwin(1:end,:).*A(1:end,:); %convolving
Z1 = fft(detrend(A));
amp_Z1 = abs(Z1).^2;
amp_Z1 = (2/(size(Z1,1)^2))*amp_Z1;
amp_Z1 = amp_Z1(1:floor(size(Z1,1)/2)+1,:);
amp_Z1 = amp_Z1/df;

%%% plotting 3 year data vs 2015 data %%%
fig = figure;

%plotting 3 year data
subplot(2,1,1)
loglog(f, nanmean(amp_Z1,2))
legend('3 year data')

%plotting 3 year and 2015 data 
subplot(2,1,2)
loglog(f, nanmean(amp_Z,2), f, nanmean(amp_Z1,2))
legend('2015 data', '3 year data')
xlabel('frequency in cpd')
h=axes(fig,'visible','off'); 
h.YLabel.Visible = 'on';
h.Title.Visible = 'on';
ylabel(h, 'Spectrum <X>^2/df in m^2/s^2 per cpd')
title('Spectrum of 3 year vs 2015 data')

%%% uncertainty for 3 vs 1 year %%%
fig = figure

%uncertainty for 3 year data
subplot(2,1,1)
nu=2*size(Z,2);                   %defining degrees of freedom
err_low = nu/chi2inv(1-.05/2,nu); %interval based on chi2inv function
err_high = nu/chi2inv(.05/2,nu);
exp = 2*sum(abs(Z(1:floor(size(Z,1)/2)+1,:)).^2,2)/size(Z,2); %calculating E(f) as sum of squares over f
exp = exp(2:end);                 %removing zero frequency
semilogy(f(2:end),exp)            %plotting E(f)
hold on 
semilogy([f(floor(length(f)/2)) f(floor(length(f)/2))], [err_low err_high]*exp(floor(length(f)/2)), 'LineWidth', 3)
xlim([0 70])
legend('2015', 'uncertainty')

%uncertainty for 2015 data
subplot(2,1,2)
nu=2*size(Z1,2);
err_low = nu/chi2inv(1-.05/2,nu);
err_high = nu/chi2inv(.05/2,nu);
exp = 2*sum(abs(Z1(1:floor(size(Z1,1)/2)+1,:)).^2,2)/size(Z1,2);
exp = exp(2:end);
semilogy(f(2:end),exp)
hold on 
semilogy([f(floor(length(f)/2)) f(floor(length(f)/2))], [err_low err_high]*exp(floor(length(f)/2)), 'LineWidth', 3)
xlim([0 70])
xlabel('frequency in cpd')
legend('3 year', 'uncertainty')

h=axes(fig,'visible','off'); 
h.YLabel.Visible = 'on';
h.Title.Visible = 'on';
ylabel(h, 'Expectation value of chi squared variable')
title('Uncertainty estimates at 95% confidence level')

%% Spectra for wspd, uwnd, vwnd

Cwspd = horzcat(Cwspd2015, Cwspd2016, Cwspd2017); %creating single wspd array
Cuwnd = horzcat(Cuwnd2015, Cuwnd2016, Cuwnd2017); %creating single uwnd array
Cvwnd = horzcat(Cvwnd2015, Cvwnd2016, Cvwnd2017); %creating single vwnd array
Ctime = horzcat(Ctime2015, Ctime2016, Ctime2017); %creating single time array

hanwin = cos(pi*Ctime/ L).^2;                     %Hanning window

A = Cwspd;
A(1:end,:) = hanwin(1:end,:).*A(1:end,:);         %convolving
Z1 = fft(detrend(A));                             %calculating spectrum as before
amp_Z1 = abs(Z1).^2;
amp_Z1 = (2/(size(Z1,1)^2))*amp_Z1;
amp_Z1 = amp_Z1(1:floor(size(Z1,1)/2)+1,:);
amp_Z1 = amp_Z1/df;

B = Cuwnd;
B(1:end,:) = hanwin(1:end,:).*B(1:end,:);
Z2 = fft(detrend(B));
amp_Z2 = abs(Z2).^2;
amp_Z2 = (2/(size(Z2,1)^2))*amp_Z2;
amp_Z2 = amp_Z2(1:floor(size(Z2,1)/2)+1,:);
amp_Z2 = amp_Z2/df;

C = Cvwnd;
C(1:end,:) = hanwin(1:end,:).*C(1:end,:);
Z3 = fft(detrend(C));
amp_Z3 = abs(Z3).^2;
amp_Z3 = (2/(size(Z3,1)^2))*amp_Z3;
amp_Z3 = amp_Z3(1:floor(size(Z3,1)/2)+1,:);
amp_Z3 = amp_Z3/df;

%%% plotting the 3 year spectra for all variables %%%
fig = figure;

subplot(3,1,1)
loglog(f,nanmean(amp_Z1,2))
legend('wind speed')
subplot(3,1,2)
loglog(f,nanmean(amp_Z2,2))
legend('u wind')
subplot(3,1,3)
loglog(f,nanmean(amp_Z3,2))
legend('v wind')
xlabel('frequency in cpd')

h=axes(fig,'visible','off'); 
h.YLabel.Visible = 'on';
h.Title.Visible = 'on';
ylabel(h, 'Spectral power in m^2/s^2 per cpd')
title('Spectra for wspd, uwnd, vwnd variables')

