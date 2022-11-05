%%%% SIOC 221A Homework 5

%% Load Data

%time [days since 1950]
TAOS15.time = ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/oceansites-tao/T8S110W/OS_T8S110W_DM134A-20150425_D_WIND_10min.nc','TIME');
TAOS16.time = ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/oceansites-tao/T8S110W/OS_T8S110W_DM183A-20160321_D_WIND_10min.nc','TIME');
TAOS17.time = ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/oceansites-tao/T8S110W/OS_T8S110W_DM231A-20170606_D_WIND_10min.nc','TIME');

%convert to datenum
TAOS15.dnum = TAOS15.time+datenum(1950,1,1);
TAOS16.dnum = TAOS16.time+datenum(1950,1,1);
TAOS17.dnum = TAOS17.time+datenum(1950,1,1);

%wind speed [m/s]
TAOS15.wspd = ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/oceansites-tao/T8S110W/OS_T8S110W_DM134A-20150425_D_WIND_10min.nc','WSPD')';
TAOS16.wspd = ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/oceansites-tao/T8S110W/OS_T8S110W_DM183A-20160321_D_WIND_10min.nc','WSPD')';
TAOS17.wspd = ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/oceansites-tao/T8S110W/OS_T8S110W_DM231A-20170606_D_WIND_10min.nc','WSPD')';

%zonal wind
TAOS15.u = ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/oceansites-tao/T8S110W/OS_T8S110W_DM134A-20150425_D_WIND_10min.nc','UWND')';
TAOS16.u = ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/oceansites-tao/T8S110W/OS_T8S110W_DM183A-20160321_D_WIND_10min.nc','UWND')';
TAOS17.u = ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/oceansites-tao/T8S110W/OS_T8S110W_DM231A-20170606_D_WIND_10min.nc','UWND')';

%meridional wind
TAOS15.v = ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/oceansites-tao/T8S110W/OS_T8S110W_DM134A-20150425_D_WIND_10min.nc','VWND')';
TAOS16.v = ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/oceansites-tao/T8S110W/OS_T8S110W_DM183A-20160321_D_WIND_10min.nc','VWND')';
TAOS17.v = ncread('https://dods.ndbc.noaa.gov/thredds/dodsC/oceansites-tao/T8S110W/OS_T8S110W_DM231A-20170606_D_WIND_10min.nc','VWND')';


%% Make a preliminary assessment of the data

%plot wind speed, zonal wind and meridional wind
figure()
plot(TAOS15.dnum,TAOS15.wspd)
hold on
plot(TAOS16.dnum,TAOS16.wspd)
plot(TAOS17.dnum,TAOS17.wspd)
title('Wind Speed')
datetick('x','mmm')
ylabel('[m/s]')
legend('2015','2016','2017')

%plot time between samples
figure()
plot(TAOS15.dnum,[mean(diff(TAOS15.time));diff(TAOS15.time)].*24.*3600)
hold on
plot(TAOS16.dnum,[mean(diff(TAOS16.time));diff(TAOS16.time)].*24.*3600)
plot(TAOS17.dnum,[mean(diff(TAOS17.time));diff(TAOS17.time)].*24.*3600)
title('Time Between Samples')
datetick('x','mmm')
ylabel('[days]')
legend('2015','2016','2017')

%plot gaps between records
figure()
plot([1,2],[TAOS16.time(end)-TAOS15.time(1),TAOS17.time(end)-TAOS16.time(1)],'o','LineWidth',2)
title('Gaps Between Records')
ylabel('Time [hours]')
xlims([0 3])

%% Linearly Interpolate to fill gaps within datsets
% There are a few sections in the data that are NaNed

% remove trailing NaNs from end of 2015 file
TAOS15.time(47582:end)=[];
TAOS15.dnum(47582:end)=[];
TAOS15.wspd(47582:end)=[];
TAOS15.u(47582:end)=[];
TAOS15.v(47582:end)=[];

%%%%%linearly interpolate over NaNs
% TAOS15
naned = isnan(TAOS15.wspd);
TAOS15.wspd(naned) = interp1(TAOS15.time(~naned),TAOS15.wspd(~naned),TAOS15.time(naned));
TAOS15.u(naned) = interp1(TAOS15.time(~naned),TAOS15.u(~naned),TAOS15.time(naned));
TAOS15.v(naned) = interp1(TAOS15.time(~naned),TAOS15.v(~naned),TAOS15.time(naned));

% TAOS16
naned = isnan(TAOS16.wspd);
TAOS16.wspd(naned) = interp1(TAOS16.time(~naned),TAOS16.wspd(~naned),TAOS16.time(naned));
TAOS16.u(naned) = interp1(TAOS16.time(~naned),TAOS16.u(~naned),TAOS16.time(naned));
TAOS16.v(naned) = interp1(TAOS16.time(~naned),TAOS16.v(~naned),TAOS16.time(naned));

% TAOS17
naned = isnan(TAOS17.wspd);
TAOS17.wspd(naned) = interp1(TAOS17.time(~naned),TAOS17.wspd(~naned),TAOS17.time(naned));
TAOS17.u(naned) = interp1(TAOS17.time(~naned),TAOS17.u(~naned),TAOS17.time(naned));
TAOS17.v(naned) = interp1(TAOS17.time(~naned),TAOS17.v(~naned),TAOS17.time(naned));



%% Segment the data into 60-day segments with 50% overlap.
%number of indices that fit into a 60 day chunk
chunk = find(TAOS15.time-TAOS15.time(1)>60,1,'first');

%TAOS15
ind1 = 1;
for i=1:10
    TAOS15.wspd1(:,i) = TAOS15.wspd(ind1:ind1+chunk);
    TAOS15.u1(:,i) = TAOS15.u(ind1:ind1+chunk);
    TAOS15.v1(:,i) = TAOS15.v(ind1:ind1+chunk);
    ind1 = ind1+chunk/2;
end
%TAOS16
ind1 = 1;
for i=1:13
    TAOS16.wspd1(:,i) = TAOS16.wspd(ind1:ind1+chunk);
    TAOS16.u1(:,i) = TAOS16.u(ind1:ind1+chunk);
    TAOS16.v1(:,i) = TAOS16.v(ind1:ind1+chunk);
    ind1 = ind1+chunk/2;
end
%TAOS17
ind1 = 1;
for i=1:9
    TAOS17.wspd1(:,i) = TAOS17.wspd(ind1:ind1+chunk);
    TAOS17.u1(:,i) = TAOS17.u(ind1:ind1+chunk);
    TAOS17.v1(:,i) = TAOS17.v(ind1:ind1+chunk);
    ind1 = ind1+chunk/2;
end


%% Determine the frequencies (in cycles per day) that you will be able to analyze when you Fourier transform the data.
dt = nanmean(diff(TAOS15.time)); %time between samples [days]
fn = 1/2/dt;  % Nyquist frequency

%TAOS15
N15 = length(TAOS15.wspd1(:,1));  %number of samples
T15 = dt*N15;   %length of record [days]
df15 = 1/T15;   % fundamental frequency [cpd]
f15 = 0:df15:fn; % frequency vector [cpd]

%TAOS16
N16 = length(TAOS16.wspd1(:,1));
T16 = dt*N16;
df16 = 1/T16;
f16 = 0:df16:fn;

%TAOS17
N17 = length(TAOS17.wspd1(:,1));
T17 = dt*N17;
df17 = 1/T17;
f17 = 0:df17:fn;

disp(strcat(['2015 frequency range: ', num2str(df15), ' cpd to ', num2str(fn), ' cpd']))
disp(strcat(['2016 frequency range: ', num2str(df16), ' cpd to ', num2str(fn), ' cpd']))
disp(strcat(['2017 frequency range: ', num2str(df17), ' cpd to ', num2str(fn), ' cpd']))

%%  Compute and plot spectra using 3 different approaches, for the 2015 wind speed record only.
% Compute the spectrum from the raw data,
for i=1:length(TAOS15.wspd1(1,:))
    a = fft(TAOS15.wspd1(:,i));
    amp=abs(a(1:(N15+1)/2)).^2; %take half of spectrum and square
    amp = amp / N15.^2; % MATLAB normalization
    amp = amp .* 2; % lost variance
    amp = amp / df15; % definition of the spectrum
    A_raw(:,i) = amp;
end
a_raw = mean(A_raw,2);

% from the detrended data,
for i=1:length(TAOS15.wspd1(1,:))
    TAOS15.wspd_dtr(:,i) = TAOS15.wspd1(:,i)-mean(TAOS15.wspd1(:,i));
    a = fft(TAOS15.wspd_dtr(:,i));
    amp=abs(a(1:(N15+1)/2)).^2; %take half of spectrum and square
    amp = amp / N15.^2; % MATLAB normalization
    amp = amp .* 2; % lost variance
    amp = amp / df15; % definition of the spectrum
    A_dtr(:,i) = amp;
end
a_dtr = mean(A_dtr,2);

% and from the detrended data with a Hanning window applied. 
%%%%% write a function for this one %%%%%
for i=1:length(TAOS15.wspd1(1,:))
    %multiply by normalized hanning window
    TAOS15.wspd_Hdtr(:,i) = TAOS15.wspd_dtr(:,i).*(hann(chunk+1).*sqrt(8/3));
    a = fft(TAOS15.wspd_Hdtr(:,i));
    amp=abs(a(1:(N15+1)/2)).^2; %take half of spectrum and square
    amp = amp / N15.^2; % MATLAB normalization
    amp = amp .* 2; % lost variance.
    amp = amp / df15; % definition of the spectrum
    A_Hdtr(:,i) = amp;
end
TAOS15.a_Hdtr = mean(A_Hdtr,2);

figure()
loglog(f15,a_raw)
hold on
loglog(f15, a_dtr)
loglog(f15, TAOS15.a_Hdtr)
title('2015 Wind Speed Spectrum')
xlabel('Frequecy [cpd]')
ylabel('\Phi_v (m/s)^2/cpd')
legend('Raw','Detrended','Detreanded and Hanning Windowed')


%% Add uncertainty estimates to the 2015 wind speed spectra.
%uncertainty with no windowing
nu_nowin=2*length(TAOS15.wspd1(1,:));
err_high15_nw = nu_nowin/chi2inv(.05/2,nu_nowin);
err_low15_nw =  nu_nowin/chi2inv(1-.05/2,nu_nowin);
%uncertainty with Hanning window
nu=0.9*2*length(TAOS15.wspd1(1,:));
err_high15 = nu/chi2inv(.05/2,nu);
err_low15 =  nu/chi2inv(1-.05/2,nu);

figure()
loglog(f15, a_raw)
hold on
loglog(f15, a_dtr)
loglog(f15, TAOS15.a_Hdtr)
loglog([0.9 0.9],[err_low15_nw err_high15_nw].*TAOS15.a_Hdtr(length(TAOS15.a_Hdtr)),'-*')
loglog([1.1 1.1],[err_low15 err_high15].*TAOS15.a_Hdtr(length(TAOS15.a_Hdtr)),'-*')
title('2015 Wind Speed Spectrum')
legend('Raw','Detrended','Detrended and Hanning Window','Uncertainty- No Window','Uncertainty- Hanning Window')
xlabel('Frequecy [cpd]')
ylabel('\Phi_v (m/s)^2/cpd')


%% Compute the wind speed spectrum using data from all three data files.
%Overlay the multi-year spectrum over the 2015 spectrum you computed earlier.

% 2016 wind speed spectrum
for i=1:length(TAOS16.wspd1(1,:))  % for each 60-day chunk
    %detrend multiply by normalized hanning window
    TAOS16.wspd_Hdtr(:,i) = (TAOS16.wspd1(:,i)-mean(TAOS16.wspd1(:,i))).*(hann(chunk+1).*sqrt(8/3));
    a = fft(TAOS16.wspd_Hdtr(:,i));  %fourier transform
    amp=abs(a(1:(N16+1)/2)).^2; %take half of spectrum and square
    amp = amp / N16.^2; % MATLAB normalization
    amp = amp .* 2; % lost variance.
    amp = amp / df16; % definition of the spectrum
    A_Hdtr(:,i) = amp;
end
TAOS16.a_Hdtr = mean(A_Hdtr,2);
% 2016 error
nu=0.9*2*length(TAOS16.wspd1(1,:));
err_high16 = nu/chi2inv(.05/2,nu);
err_low16 =  nu/chi2inv(1-.05/2,nu);


% 2017 wind speed spectrum
for i=1:length(TAOS17.wspd1(1,:))  % for each 60-day chunk
    %detrend multiply by normalized hanning window
    TAOS17.wspd_Hdtr(:,i) = (TAOS17.wspd1(:,i)-mean(TAOS17.wspd1(:,i))).*(hann(chunk+1).*sqrt(8/3));
    a = fft(TAOS17.wspd_Hdtr(:,i));  %fourier transform
    amp=abs(a(1:(N17+1)/2)).^2; %take half of spectrum and square
    amp = amp / N17.^2; % MATLAB normalization
    amp = amp .* 2; % lost variance.
    amp = amp / df17; % definition of the spectrum
    A_Hdtr(:,i) = amp;
end
TAOS17.a_Hdtr = mean(A_Hdtr,2);
% 2017 uncertainty
nu=0.9*2*length(TAOS17.wspd1(1,:));
err_high17 = nu/chi2inv(.05/2,nu);
err_low17 =  nu/chi2inv(1-.05/2,nu);

figure()
loglog(f15,TAOS15.a_Hdtr,'-r')
hold on
loglog([0.8 0.8],[err_low15 err_high15].*4e-3,'-*r')
loglog(f16,TAOS16.a_Hdtr,'-b')
loglog([1 1],[err_low16 err_high16].*4e-3,'-*b')
loglog(f17,TAOS17.a_Hdtr,'-g')
loglog([1.2 1.2],[err_low17 err_high17].*4e-3,'-*g')
loglog([1 1],[min(TAOS15.a_Hdtr) max(TAOS15.a_Hdtr)],'--k')
loglog([2 2],[min(TAOS15.a_Hdtr) max(TAOS15.a_Hdtr)],'--k')
title('Wind Speed Spectra')
legend('2015','-','2016','-','2017','-')
xlabel('Frequecy [cpd]')
ylabel('\Phi_v (m/s)^2/cpd')


%% Compare the spectra for wind speed, zonal wind, and meridional wind.


%%%%%%%%%%%%%%%%%%%%%%%%%%%2015
%zonal wind spectrum
[~,TAOS15.au_Hdtr] = spectrumCB(TAOS15.time,TAOS15.u,chunk,1,1);
% for i=1:length(TAOS15.u1(1,:))  % for each 60-day chunk
%     %detrend multiply by normalized hanning window
%     TAOS15.u_Hdtr(:,i) = (TAOS15.u1(:,i)-mean(TAOS15.u1(:,i))).*(hann(chunk+1)./mean(hann(chunk+1)));
%     a = fft(TAOS15.u_Hdtr(:,i));  %fourier transform
%     amp=abs(a(1:(N15+1)/2)).^2; %take half of spectrum and square
%     amp = amp / N15.^2; % MATLAB normalization
%     amp = amp .* 2; % lost variance.
%     amp = amp / df15; % definition of the spectrum
%     A_Hdtr(:,i) = amp;
% end
% TAOS15.au_Hdtr = mean(A_Hdtr,2);

%meridional wind spectrum
[f15,TAOS15.av_Hdtr] = spectrumCB(TAOS15.time,TAOS15.v,chunk,1,1);
% for i=1:length(TAOS15.v1(1,:))  % for each 60-day chunk
%     %detrend multiply by normalized hanning window
%     TAOS15.v_Hdtr(:,i) = (TAOS15.v1(:,i)-mean(TAOS15.v1(:,i))).*(hann(chunk+1)./mean(hann(chunk+1)));
%     a = fft(TAOS15.v_Hdtr(:,i));  %fourier transform
%     amp=abs(a(1:(N15+1)/2)).^2; %take half of spectrum and square
%     amp = amp / N15.^2; % MATLAB normalization
%     amp = amp .* 2; % lost variance.
%     amp = amp / df15; % definition of the spectrum
%     A_Hdtr(:,i) = amp;
% end
% TAOS15.av_Hdtr = mean(A_Hdtr,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2016
%zonal wind spectrum
[~,TAOS16.au_Hdtr] = spectrumCB(TAOS16.time,TAOS16.u,chunk,1,1);
% for i=1:length(TAOS16.u1(1,:))  % for each 60-day chunk
%     %detrend multiply by normalized hanning window
%     TAOS16.u_Hdtr(:,i) = (TAOS16.u1(:,i)-mean(TAOS16.u1(:,i))).*(hann(chunk+1)./mean(hann(chunk+1)));
%     a = fft(TAOS16.u_Hdtr(:,i));  %fourier transform
%     amp=abs(a(1:(N16+1)/2)).^2; %take half of spectrum and square
%     amp = amp / N16.^2; % MATLAB normalization
%     amp = amp .* 2; % lost variance.
%     amp = amp / df16; % definition of the spectrum
%     A_Hdtr(:,i) = amp;
% end
% TAOS16.au_Hdtr = mean(A_Hdtr,2);

%meridional wind spectrum
[f16,TAOS16.av_Hdtr] = spectrumCB(TAOS16.time,TAOS16.v,chunk,1,1);
% for i=1:length(TAOS16.v1(1,:))  % for each 60-day chunk
%     %detrend multiply by normalized hanning window
%     TAOS16.v_Hdtr(:,i) = (TAOS16.v1(:,i)-mean(TAOS16.v1(:,i))).*(hann(chunk+1)./mean(hann(chunk+1)));
%     a = fft(TAOS16.v_Hdtr(:,i));  %fourier transform
%     amp=abs(a(1:(N16+1)/2)).^2; %take half of spectrum and square
%     amp = amp / N16.^2; % MATLAB normalization
%     amp = amp .* 2; % lost variance.
%     amp = amp / df16; % definition of the spectrum
%     A_Hdtr(:,i) = amp;
% end
% TAOS16.av_Hdtr = mean(A_Hdtr,2);


%%%%%%%%%%%%%%%%%%%2017
%zonal wind spectrum
[~,TAOS17.au_Hdtr] = spectrumCB(TAOS17.time,TAOS17.u,chunk,1,1);
% for i=1:length(TAOS17.u1(1,:))  % for each 60-day chunk
%     %detrend multiply by normalized hanning window
%     TAOS17.u_Hdtr(:,i) = (TAOS17.u1(:,i)-mean(TAOS17.u1(:,i))).*(hann(chunk+1)./mean(hann(chunk+1)));
%     a = fft(TAOS17.u_Hdtr(:,i));  %fourier transform
%     amp=abs(a(1:(N17+1)/2)).^2; %take half of spectrum and square
%     amp = amp / N17.^2; % MATLAB normalization
%     amp = amp .* 2; % lost variance.
%     amp = amp / df17; % definition of the spectrum
%     A_Hdtr(:,i) = amp;
% end
% TAOS17.au_Hdtr = mean(A_Hdtr,2);

%meridional wind spectrum
[f17,TAOS17.av_Hdtr] = spectrumCB(TAOS17.time,TAOS17.v,chunk,1,1);
% for i=1:length(TAOS17.v1(1,:))  % for each 60-day chunk
%     %detrend multiply by normalized hanning window
%     TAOS17.v_Hdtr(:,i) = (TAOS17.v1(:,i)-mean(TAOS17.v1(:,i))).*(hann(chunk+1)./mean(hann(chunk+1)));
%     a = fft(TAOS17.v_Hdtr(:,i));  %fourier transform
%     amp=abs(a(1:(N17+1)/2)).^2; %take half of spectrum and square
%     amp = amp / N17.^2; % MATLAB normalization
%     amp = amp .* 2; % lost variance.
%     amp = amp / df17; % definition of the spectrum
%     A_Hdtr(:,i) = amp;
% end
% TAOS17.av_Hdtr = mean(A_Hdtr,2);

%%%% 2015 plotting
figure()
loglog(f15,TAOS15.a_Hdtr)
hold on
loglog(f15,TAOS15.au_Hdtr)
loglog(f15,TAOS15.av_Hdtr)
loglog([1 1],[min(TAOS15.a_Hdtr) max(TAOS15.a_Hdtr)],'--')
loglog([2 2],[min(TAOS15.a_Hdtr) max(TAOS15.a_Hdtr)],'--')
loglog([0.8 0.8],[err_low15 err_high15].*4e-3,'-*k')
title('2015 Wind Spectra')
legend('Wind Speed','Zonal Wind','Meridional Wind','Diurnal','Semidiurnal','Uncertainty')
xlabel('Frequecy [cpd]')
ylabel('\Phi_v (m/s)^2/cpd')

%%% 2016 plotting
figure()
loglog(f16,TAOS16.a_Hdtr)
hold on
loglog(f16,TAOS16.au_Hdtr)
loglog(f16,TAOS16.av_Hdtr)
loglog([1 1],[min(TAOS16.a_Hdtr) max(TAOS16.a_Hdtr)],'--')
loglog([2 2],[min(TAOS16.a_Hdtr) max(TAOS16.a_Hdtr)],'--')
loglog([0.8 0.8],[err_low16 err_high16].*4e-3,'-*k')
title('2016 Wind Spectra')
legend('Wind Speed','Zonal Wind','Meridional Wind','Diurnal','Semidiurnal','Uncertainty')
xlabel('Frequecy [cpd]')
ylabel('\Phi_v (m/s)^2/cpd')


%%% 2017 plotting
figure()
loglog(f17,TAOS17.a_Hdtr)
hold on
loglog(f17,TAOS17.au_Hdtr)
loglog(f17,TAOS17.av_Hdtr)
loglog([1 1],[min(TAOS17.a_Hdtr) max(TAOS17.a_Hdtr)],'--')
loglog([2 2],[min(TAOS17.a_Hdtr) max(TAOS17.a_Hdtr)],'--')
loglog([0.8 0.8],[err_low17 err_high17].*4e-3,'-*k')
title('2017 Wind Spectra')
legend('Wind Speed','Zonal Wind','Meridional Wind','Diurnal','Semidiurnal','Uncertainty')
xlabel('Frequecy [cpd]')
ylabel('\Phi_v (m/s)^2/cpd')






