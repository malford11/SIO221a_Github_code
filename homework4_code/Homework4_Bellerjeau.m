%%% SIOC 221A Homework 3

%% Load Data
%load 2021 pressure and time
time0 = ncread('PierData\scripps_pier-2021.nc','time');
pressure0 = ncread('PierData\scripps_pier-2021.nc','pressure');

%editing gaps in data, thank you Matthew
startline=70521;
endline=88190;
xx=startline:endline;
N=(time0(endline)-time0(startline))/240;
time = double(time0(startline))+240*(0:N);
time = double(time);
dnum = double(time)/86400+datenum(1970,1,1);
% gaps at 2420, 4952, 7496, 15042
pressure(1:2420)=pressure0(xx(1:2420));
pressure(2421)=0.5*sum(pressure0(xx(2420:2421)));
pressure(2422:4953)=pressure0(xx(2421:4952));
pressure(4954)=0.5*sum(pressure0(xx(4952:4953)));
pressure(4955:7498)=pressure0(xx(4953:7496));
pressure(7499)=0.5*sum(pressure0(xx(7496:7497)));
pressure(7500:15045)=pressure0(xx(7497:15042));
pressure(15046)=0.5*sum(pressure0(xx(15042:15043)));
pressure(15047:length(xx)+4)=pressure0(xx(15043:end));
pressure=double(pressure);


%% Visual Evaluation
figure()
plot(dnum, pressure)
title('2021 Pressure')
datetick('x','mmm')
xlims([dnum(1),dnum(end)])
ylabel('Days')

figure()
plot(dnum,[mean(diff(time)), diff(time)])
title('Time Between Samples')
ylabel('Seconds')
datetick('x','mmm')
xlims([dnum(1),dnum(end)])

figure()
plot(double(time0)/86400+datenum(1970,1,1),[mean(diff(time0)); diff(time0)])
title('Whole Year 2021 Time Between Samples')
ylabel('Seconds')
datetick('x','mmm')
ylims([-10 2000])



%% Least Squares Fit
f = [1/25.82*24,1/23.93*24,1/12.42*24];

[pA,px] = LSfit(time'/86400,pressure',f);

figure()
plot(dnum,pressure,dnum,pA*px)
datetick('x','mmm')
xlims([dnum(1),dnum(end)])
ylabel('Pressure [dbar]')
legend('Raw Data','Model')
title('2021 Pressure')

%amplitudes of tidal constituents
mn = px(1);
amp1 = sqrt(px(2)^2+px(3)^2);
amp2 = sqrt(px(4)^2+px(5)^2);
amp3 = sqrt(px(6)^2+px(7)^2);

%display tidal amplitudes
disp(strcat(['Mean pressure: ', num2str(mn)]))
disp(strcat(['Tidal Amplitude 1: ', num2str(amp1)]))
disp(strcat(['Tidal Amplitude 2: ', num2str(amp2)]))
disp(strcat(['Tidal Amplitude 3: ', num2str(amp3)]))


%% Fourier Transform

Fp = fftshift(fft(pressure));
dt = 0.0028;
fs = 1/dt; % sampling frequency in 1/days
T = length(pressure)*dt; %length of record in days
df = 1/T; % spacing for frequency vecotr
freqFp = -fs/2:df:fs/2; % generate frequency vector

figure()
subplot(1,2,1)
plot(freqFp(1:end-1),abs(Fp))
xlabel('Frequency [cpd]')
ylabel('|fft(pressure)|')
xlims([-4 4])
title('Fourier Transform of Pressure')
subplot(1,2,2)
plot(freqFp(1:end-1),abs(Fp))
xlabel('Frequency [cpd]')
ylabel('|fft(pressure)|')
xlims([-4 4])
ylims([0 5000])
title('Zoomed Fourier Transform of Pressure')

% amplitudes
mn_f = 60301.9;
amp1_f = 1685.98*2;
amp2_f = 2996.2*2;
amp3_f = 4134.48*2;

figure()
subplot(1,2,1)
plot(freqFp(1:end-1),abs(Fp)./mn_f)
xlabel('Frequency [cpd]')
ylabel('|fft(pressure)|./amplitude of mean')
xlims([-4 4])
title('Scaled Fourier Transform of Pressure')
subplot(1,2,2)
plot(freqFp(1:end-1),abs(Fp)./mn_f)
xlabel('Frequency [cpd]')
ylabel('|fft(pressure)|./amplitude of mean')
xlims([-4 4])
ylims([0 0.1])
title('Zoomed Scaled Fourier Transform of Pressure')



% spectral energy
Fpp = fftshift(fft(pressure-mean(pressure)));

figure()
plot(freqFp(1:end-1),(abs(Fpp).^2)./df)
xlims([-4 4])
title('Spectral Energy Density')
ylabel('|fft(pressure)|^2/df')
xlabel('Frequency [cpd]')





