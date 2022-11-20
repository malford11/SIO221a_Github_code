% Homework 6 
% SIOC 221A
% By Charlotte Bellerjeau 11/8/22

%% Compute 2 spectra

%generate random data
W = randn(10000,1);   %white gaussian noise
dt = 1; %sampling rate [1 second]
R = cumsum(W)*dt;   %integral of white noise
N = length(W);   %length of record
time = 1:dt:N;    %time vector [seconds]

%plot random data
figure()
subplot(2,1,1)
plot(time,W)
title('White Noise Data')
xlabel('Time [s]')
ylabel('Data [units]')
subplot(2,1,2)
plot(time,R)
title('Integral of White Noise Data')
xlabel('Time [s]')
ylabel('Data [units]')

%compute and plot spectra with 50% overlap and no windowing
chunk=500;  %chunk length
detrend = 1;
Hanning = 0;
cosine = 0;
[Wf,Wspec,WParseval] = spectrumCB(time,W,chunk,detrend,Hanning,cosine);
[Rf,Rspec,RParseval] = spectrumCB(time,R,chunk,detrend,Hanning,cosine);
figure()
subplot(2,1,1)
loglog(Wf,Wspec)
title('White Noise Spectrum')
xlabel('Frequency [cps]')
ylabel('\Phi_v (m/s)^2/cps')
subplot(2,1,2)
loglog(Rf,Rspec)
title('Spectrum of Integral of White Noise')
xlabel('Frequency [cps]')
ylabel('\Phi_v (m/s)^2/cps')

%% Analytical Spectra
%overplot the expected analytical spectrum on the white noise spectrum
Wvar= mean(W.^2);  %variance
fn = 1/2/dt;  % Nyquist frequency
Wlevel = Wvar./fn;
Waspec = Wlevel.*ones(length(Wspec),1);
Raspec = Waspec./(2*pi*1i.*Rf);

figure()
subplot(2,1,1)
loglog(Wf,Wspec)
hold on
plot(Wf,Waspec)
title('White Noise Spectrum')
xlabel('Frequency [cps]')
ylabel('\Phi_v (m/s)^2/cps')
legend('Spectrum','Analytical Spectrum')
subplot(2,1,2)
loglog(Rf,Rspec)
hold on
plot(Rf,abs(Raspec))
title('Spectrum of Integral of White Noise')
xlabel('Frequency [cps]')
ylabel('\Phi_v (m/s)^2/cps')

%% Error Bars
%uncertainty with no windowing
nu=2*38;
err_high = nu/chi2inv(.05/2,nu);
err_low =  nu/chi2inv(1-.05/2,nu);

figure()
subplot(2,1,1)
loglog(Wf,Wspec)
hold on
plot(Wf,Waspec)
loglog([2e-2 2e-2],[err_low err_high].*1.9827,'-*')
legend('Spectrum','Analytical Spectrum','Error Bars')
title('White Noise Spectrum')
xlabel('Frequency [cps]')
ylabel('\Phi_v (m/s)^2/cps')
subplot(2,1,2)
loglog(Rf,Rspec)
hold on
plot(Rf,abs(Raspec))
loglog([2e-2 2e-2],[err_low err_high].*1.9827,'-*')
title('Spectrum of Integral of White Noise')
xlabel('Frequency [cps]')
ylabel('\Phi_v (m/s)^2/cps')


%% Use Monte Carlo simulation to verify the χ^2 error bars
Hanning = 0;
cosine = 0;
[err_lowMC,err_highMC] = MonteCarlo_errbars(Hanning,cosine);

MCrat = mean(err_lowMC)/mean(err_highMC);
chi2rat = err_low/err_high;
disp(strcat(['The ratio of the Monte Carlo error bars is: ',num2str(MCrat)]))
disp(strcat(['The ratio of the Chi^2 error bars was: ',num2str(chi2rat)]))

figure()
loglog(Wf,Wspec)
hold on
plot(Wf,Waspec)
loglog([2e-2 2e-2],[err_low err_high].*1.9827,'-*')
plot(Wf,err_lowMC,'--','Color',[0.5 0.5 0.5])
plot(Wf,err_highMC,'--','Color',[0.5 0.5 0.5])
legend('Spectrum','Analytical Spectrum','Error Bars','Monte Carlo Error Bars','-')
title('White Noise Spectrum')
xlabel('Frequency [cps]')
ylabel('\Phi_v (m/s)^2/cps')



%% Monte Carlo process for a Hanning window
Hanning = 1;
cosine = 0;
[err_lowMC,err_highMC] = MonteCarlo_errbars(Hanning,cosine);

MCrat = mean(err_lowMC)/mean(err_highMC);
chi2rat = err_low/err_high;
disp(strcat(['The ratio of the Monte Carlo error bars is: ',num2str(MCrat)]))
disp(strcat(['The ratio of the Chi^2 error bars was: ',num2str(chi2rat)]))

figure()
loglog(Wf,Wspec)
hold on
plot(Wf,Waspec)
loglog([2e-2 2e-2],[err_low err_high].*1.9827,'-*')
plot(Wf,err_lowMC,'--','Color',[0.5 0.5 0.5])
plot(Wf,err_highMC,'--','Color',[0.5 0.5 0.5])
legend('Spectrum','Analytical Spectrum','Error Bars','Monte Carlo Error Bars','-')
title('White Noise Spectrum - Hanning Windowed')
xlabel('Frequency [cps]')
ylabel('\Phi_v (m/s)^2/cps')

%% Monte Carlo process for a Cosine window. 
Hanning = 0;
cosine = 1;
[err_lowMC,err_highMC] = MonteCarlo_errbars(Hanning,cosine);

MCrat = mean(err_lowMC)/mean(err_highMC);
chi2rat = err_low/err_high;
disp(strcat(['The ratio of the Monte Carlo error bars is: ',num2str(MCrat)]))
disp(strcat(['The ratio of the Chi^2 error bars was: ',num2str(chi2rat)]))

figure()
loglog(Wf,Wspec)
hold on
plot(Wf,Waspec)
loglog([2e-2 2e-2],[err_low err_high].*1.9827,'-*')
plot(Wf,err_lowMC,'--','Color',[0.5 0.5 0.5])
plot(Wf,err_highMC,'--','Color',[0.5 0.5 0.5])
legend('Spectrum','Analytical Spectrum','Error Bars','Monte Carlo Error Bars','-')
title('White Noise Spectrum - Cosine Windowed')
xlabel('Frequency [cps]')
ylabel('\Phi_v (m/s)^2/cps')
