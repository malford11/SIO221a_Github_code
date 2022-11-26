%%%%SIOC 221A Data Analysis
%%%%Homework 8
%%%%Charlotte Bellerjeau

%% Revisiting normalization and Parseval.

%%%%%%choose any real dataset from this course
%time [days since 1950]
TAOS15.time = ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/oceansites-tao/T8S110W/OS_T8S110W_DM134A-20150425_D_WIND_10min.nc','TIME');
%convert to datenum
TAOS15.dnum = TAOS15.time+datenum(1950,1,1);
%windspeed
TAOS15.wspd = ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/oceansites-tao/T8S110W/OS_T8S110W_DM134A-20150425_D_WIND_10min.nc','WSPD')';
% remove trailing NaNs from end of record
TAOS15.time(47582:end)=[];
TAOS15.dnum(47582:end)=[];
TAOS15.wspd(47582:end)=[];
%linearly interpolate over NaNs
% TAOS15
naned = isnan(TAOS15.wspd);
TAOS15.wspd(naned) = interp1(TAOS15.time(~naned),TAOS15.wspd(~naned),TAOS15.time(naned));



%%%%%%show that parseval's theorem is satisfied for 2 degrees of freedom

% frequency vector
dt = nanmean(diff(TAOS15.time)); %time between samples [days]
fn = 1/2/dt;  % Nyquist frequency
N= length(TAOS15.wspd);
T = dt*N;   %length of record [days]
df = 1/T;   % fundamental frequency [cpd]
f = 0:df:fn; % frequency vector [cpd]
f = f';
%spectrum
data = detrend(TAOS15.wspd);  %just detrend
a = fft(data);
amp=abs(a(1:(N+1)/2)).^2; %take half of spectrum and square
amp = amp / N.^2; % MATLAB normalization
amp = amp .* 2; % lost variance
amp = amp / df; % definition of the spectrum
%check Parseval's theorem
variance=std(data)^2;
int_spec = trapz(f,amp);
Pars2=int_spec/variance;

%plot spectrum
figure()
loglog(f,amp)
title('TAOS 2015 Wind Speed Spectrum, 2 DOF')
xlabel('Frequecy [cpd]')
ylabel('\Phi_v (m/s)^2/cpd')

%%%%%%show parseval for 10-20 (ok,15?) half overlapped segments with a Hanning window
chunk=6000;
dtrend=1;
Hanning=1;
cosine=0;
[f15,a15,Pars15]=spectrumCB(TAOS15.time, TAOS15.wspd, chunk, dtrend, Hanning, cosine);

%plot spectrum
figure()
loglog(f15,a15)
title('TAOS 2015 Wind Speed Spectrum, 2 DOF')
xlabel('Frequecy [cpd]')
ylabel('\Phi_v (m/s)^2/cpd')

%% Aliasing.

%sampling frequencies [cpd]
fs_fst =  1/0.99349;
fs_sci = 1/20.86455;
%tidal frequencies
fS1 = 1/24.00;
f2N2 = 1/12.9054;
fN2 = 1/12.6583;
fM2 = 1/12.4206;
fS2 = 1/12.00;
fK2 = 1/11.9672;

%%% aliasing of each tidal frequency by each orbit
%fast orbit
faS1fst = alias(fS1,fs_fst);
fa2N2fst = alias(f2N2,fs_fst);
faN2fst = alias(fN2,fs_fst);
faM2fst = alias(fM2,fs_fst);
faS2fst = alias(fS2,fs_fst);
faK2fst = alias(fK2,fs_fst);
%science orbit
faS1sci = alias(fS1,fs_sci);
fa2N2sci = alias(f2N2,fs_sci);
faN2sci = alias(fN2,fs_sci);
faM2sci = alias(fM2,fs_sci);
faS2sci = alias(fS2,fs_sci);
faK2sci = alias(fK2,fs_sci);

%%%how long should the record be to resolve all of these peaks?
dfs = fS2-fK2; %frequency difference of closest peaks (S2 and K2)
T = abs(1/dfs);


disp('The tides will be aliased by the fast orbit as below:')
disp(strcat(['S1:  ',num2str(faS1fst)]))
disp(strcat(['2N2:  ',num2str(fa2N2fst)]))
disp(strcat(['N2:  ',num2str(faN2fst)]))
disp(strcat(['M2:  ',num2str(faM2fst)]))
disp(strcat(['S2:  ',num2str(faS2fst)]))
disp(strcat(['K2:  ',num2str(faK2fst)]))

disp('The tides will be aliased by the science orbit as below:')
disp(strcat(['S1:  ',num2str(faS1sci)]))
disp(strcat(['2N2:  ',num2str(fa2N2sci)]))
disp(strcat(['N2:  ',num2str(faN2sci)]))
disp(strcat(['M2:  ',num2str(faM2sci)]))
disp(strcat(['S2:  ',num2str(faS2sci)]))
disp(strcat(['K2:  ',num2str(faK2sci)]))


disp('The closest two tidal peaks belong to the S2 and K2 tides.')
disp(strcat(['To resolve these two peaks the record must be at least ', num2str(T),' days long, regardless of the orbital period.']))


%% Spectra of aliased signals.

%compute the spectrum of the TAOS windspeed data with 60 day overlapping chunks
chunk = 7800; %53.8 day window, fits nicer in time series with 50% overlap
dtrend=1;
Hanning=1;
cosine=0;
[f,a,Pars]=spectrumCB(TAOS15.time, TAOS15.wspd, chunk, dtrend, Hanning, cosine);

%subsample at every 40 data points and recompute the spectrum
sub40 = 1:40:47560;
TAOS15.wspd_sub = TAOS15.wspd(sub40);
TAOS15.time_sub = TAOS15.time(sub40);
chunk = 200;
[fsub,asub,ParsSub]=spectrumCB(TAOS15.time_sub, TAOS15.wspd_sub, chunk, dtrend, Hanning, cosine);

figure()
loglog(f,a)
hold on
loglog(fsub,asub)
loglog([1 1],[min(a),max(a)],'--k')
loglog([2 2],[min(a),max(a)],'--k')
title('TAOS 2015 Wind Speed Spectrum')
legend('original','subsampled')
xlabel('Frequecy [cpd]')
ylabel('\Phi_v (m/s)^2/cpd')

%resolution and Nyquist frequency of subsampled data
dt = nanmean(diff(TAOS15.time_sub)); %time between samples [days]
fn = 1/2/dt;  % Nyquist frequency

disp(strcat(['The resolution of the subsampled time series is: ', num2str(dt)]))
disp(strcat(['The Nyquist frequency of the subsampled time series is: ', num2str(fn)]))

%alias frequency of semi-diurnal peak for subsampled time series
fA = alias(1/dt,2);

%how do the spectral energies compare?
E0 = Pars*std(TAOS15.wspd)^2;
Esub = ParsSub*std(TAOS15.wspd_sub)^2;


