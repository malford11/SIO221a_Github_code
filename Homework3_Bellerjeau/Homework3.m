%%%% SIOC 221A Data Analysis
%%%% Homework 3
%%%% Charlotte Bellerjeau

clear all

%% Visual Evaluation
% Load bouy data from 2015 and 2016
delim = ' ';
header = 2;
bouy2015 = importdata('bouy46047_2015.txt',delim,header);
bouy2016 = importdata('bouy46047_2016.txt',delim,header);

for i=1:length(bouy2015.data(:,1))
    bouy2015.dnum(i) = datenum([bouy2015.data(i,1),bouy2015.data(i,2),bouy2015.data(i,3),bouy2015.data(i,4),bouy2015.data(i,5),0]);
end
for i=1:length(bouy2016.data(:,1))
    bouy2016.dnum(i) = datenum([bouy2016.data(i,1),bouy2016.data(i,2),bouy2016.data(i,3),bouy2016.data(i,4),bouy2016.data(i,5),0]);
end

bouy.dnum = [bouy2015.dnum';bouy2016.dnum'];
bouy.wspd = [bouy2015.data(:,7);bouy2016.data(:,7)];
bouy.wh = [bouy2015.data(:,9);bouy2016.data(:,9)];
bouy.wtemp = [bouy2015.data(:,15);bouy2016.data(:,15)];
bouy.atemp = [bouy2015.data(:,14);bouy2016.data(:,14)];


% Plot wind speed, wave height, water temperature, and air temperature
figure()
subplot(4,1,1)
plot(bouy.dnum,bouy.wspd)
datetick('x','mmm')
title('Wind Speed')
ylabel('[m/s]')
subplot(4,1,2)
plot(bouy.dnum,bouy.wh)
datetick('x','mmm')
title('Wave Height')
ylabel('[m]')
subplot(4,1,3)
plot(bouy.dnum,bouy.wtemp)
datetick('x','mmm')
title('Water Temperature')
ylabel('[degC]')
subplot(4,1,4)
plot(bouy.dnum,bouy.atemp)
datetick('x','mmm')
title('Air Temperature')
ylabel('[degC]')

%Removing outliers and replotting
bouy.wspd(bouy.wspd>50) = bouy.wspd(find(bouy.wspd>50)-1);
bouy.wh(bouy.wh>50) = bouy.wh(find(bouy.wh>50)-1);
bouy.wh(bouy.wh>50) = bouy.wh(find(bouy.wh>50)-1);
bouy.wh(bouy.wh>50) = bouy.wh(find(bouy.wh>50)-1);
bouy.wh(bouy.wh>50) = bouy.wh(find(bouy.wh>50)-1);
bouy.wtemp(bouy.wtemp>500) = bouy.wtemp(find(bouy.wtemp>500)-1);
bouy.atemp(bouy.atemp>500) = bouy.atemp(find(bouy.atemp>500)-1);

% Plot wind speed, wave height, water temperature, and air temperature
figure()
subplot(4,1,1)
plot(bouy.dnum,bouy.wspd)
datetick('x','mmm')
title('Wind Speed')
ylabel('[m/s]')
subplot(4,1,2)
plot(bouy.dnum,bouy.wh)
datetick('x','mmm')
title('Wave Height')
ylabel('[m]')
subplot(4,1,3)
plot(bouy.dnum,bouy.wtemp)
datetick('x','mmm')
title('Water Temperature')
ylabel('[degC]')
subplot(4,1,4)
plot(bouy.dnum,bouy.atemp)
datetick('x','mmm')
title('Air Temperature')
ylabel('[degC]')


%% Monthly means
bouy.yday=datenum2yday(bouy.dnum);
bouy.yday=bouy.yday-bouy.yday(1);

months=0:30:24*30;
mind=1;
for i=1:length(months)
    mind(i+1) = find(bouy.yday>months(i),1,'first');
    bouy.wspdAvg(i) = mean(bouy.wspd(mind(i):mind(i+1)));
    bouy.whAvg(i) = mean(bouy.wh(mind(i):mind(i+1)));
    bouy.wtempAvg(i) = mean(bouy.wtemp(mind(i):mind(i+1)));
    bouy.atempAvg(i) = mean(bouy.atemp(mind(i):mind(i+1)));
    bouy.dnumAvg(i) = mean(bouy.dnum(mind(i):mind(i+1)));
    bouy.wspdAerr(i) = std(bouy.wspd(mind(i):mind(i+1)))/sqrt(30/7);
    bouy.whAerr(i) = std(bouy.wh(mind(i):mind(i+1)))/sqrt(30/7);
    bouy.wtempAerr(i) = std(bouy.wtemp(mind(i):mind(i+1)))/sqrt(30/7);
    bouy.atempAerr(i) = std(bouy.atemp(mind(i):mind(i+1)))/sqrt(30/7);
end

%Plot monthly means with errorbars
figure()
subplot(4,1,1)
errorbar(bouy.dnumAvg,bouy.wspdAvg,bouy.wspdAerr)
datetick('x','mmm')
title('Average Wind Speed')
ylabel('[m/s]')
subplot(4,1,2)
errorbar(bouy.dnumAvg,bouy.whAvg,bouy.whAerr)
datetick('x','mmm')
title('Average Wave Height')
ylabel('[m]')
subplot(4,1,3)
errorbar(bouy.dnumAvg,bouy.wtempAvg,bouy.wtempAerr)
datetick('x','mmm')
title('Average Water Temperature')
ylabel('[degC]')
subplot(4,1,4)
errorbar(bouy.dnumAvg,bouy.atempAvg,bouy.atempAerr)
datetick('x','mmm')
title('Average Air Temperature')
ylabel('[degC]')



%% Least Squares Fit
%split the data back into two seperate years for fitting
by15.dnum=bouy.dnumAvg(1:12);
by15.wspd=bouy.wspdAvg(1:12);
by15.wh=bouy.whAvg(1:12);
by15.wtemp=bouy.wtempAvg(1:12);
by15.atemp=bouy.atempAvg(1:12);
by15.wspdErr = bouy.wspdAerr(1:12);
by15.whErr = bouy.whAerr(1:12);
by15.wtempErr = bouy.wtempAerr(1:12);
by15.atempErr = bouy.atempAerr(1:12);
by15.wspdErr(1) = by15.wspdErr(2);
by15.whErr(1) = by15.wspdErr(2);
by15.wtempErr(1) = by15.wspdErr(2);
by15.atempErr(1) = by15.wspdErr(2);


by16.dnum=bouy.dnumAvg(12:end);
by16.wspd=bouy.wspdAvg(12:end);
by16.wh=bouy.whAvg(12:end);
by16.wtemp=bouy.wtempAvg(12:end);
by16.atemp=bouy.atempAvg(12:end);
by16.wspdErr = bouy.wspdAerr(12:end);
by16.whErr = bouy.whAerr(12:end);
by16.wtempErr = bouy.wtempAerr(12:end);
by16.atempErr = bouy.atempAerr(12:end);

%least squares fit the four datasets for each year with a mean and annual cycle
f=1/365;
%2015
[by15.wspdA,by15.wspdx]=LSfit(by15.dnum,by15.wspd,f);
[by15.whA,by15.whx]=LSfit(by15.dnum,by15.wh,f);
[by15.wtempA,by15.wtempx]=LSfit(by15.dnum,by15.wtemp,f);
[by15.atempA,by15.atempx]=LSfit(by15.dnum,by15.atemp,f);
%2016
[by16.wspdA,by16.wspdx]=LSfit(by16.dnum,by16.wspd,f);
[by16.whA,by16.whx]=LSfit(by16.dnum,by16.wh,f);
[by16.wtempA,by16.wtempx]=LSfit(by16.dnum,by16.wtemp,f);
[by16.atempA,by16.atempx]=LSfit(by16.dnum,by16.atemp,f);

%plot data along with LSfits
figure()
subplot(4,2,1)
plot(by15.dnum,by15.wspd,by15.dnum,by15.wspdA*by15.wspdx)
datetick('x','mmm')
title('2015 Annual Fitted Average Wind Speed')
ylabel('[m/s]')
subplot(4,2,3)
plot(by15.dnum,by15.wh,by15.dnum,by15.whA*by15.whx)
datetick('x','mmm')
title(' "" Wave Height')
ylabel('[m]')
subplot(4,2,5)
plot(by15.dnum,by15.wtemp,by15.dnum,by15.wtempA*by15.wtempx)
datetick('x','mmm')
title(' "" Water Temperature')
ylabel('[degC]')
subplot(4,2,7)
plot(by15.dnum,by15.atemp,by15.dnum,by15.atempA*by15.atempx)
datetick('x','mmm')
title(' "" Air Temperature')
ylabel('[degC]')


subplot(4,2,2)
plot(by16.dnum,by16.wspd,by16.dnum,by16.wspdA*by16.wspdx)
datetick('x','mmm')
title('2016 Annual Fitted Average Wind Speed')
ylabel('[m/s]')
subplot(4,2,4)
plot(by16.dnum,by16.wh,by16.dnum,by16.whA*by16.whx)
datetick('x','mmm')
title(' "" Wave Height')
ylabel('[m]')
subplot(4,2,6)
plot(by16.dnum,by16.wtemp,by16.dnum,by16.wtempA*by16.wtempx)
datetick('x','mmm')
title(' "" Water Temperature')
ylabel('[degC]')
subplot(4,2,8)
plot(by16.dnum,by16.atemp,by16.dnum,by16.atempA*by16.atempx)
datetick('x','mmm')
title(' "" Air Temperature')
ylabel('[degC]')

% display mean and amplitude of the annual cycle for each year
disp('______________')
disp('2015')
disp('______________')
disp('--Annual Fit Means--')
disp('Wind Speed: 6.05 m/s')
disp('Wave Height: 1.98 m')
disp('Water Temperature: 17.99 C')
disp('Air Temperature: 16.84 C')
disp('--Annual Cycle Amplitude--')
disp(strcat(['Wind Speed: ', num2str(sqrt(sum(by15.wspdx(2:end).^2))), '  m/s']))
disp(strcat(['Wave Height: ', num2str(sqrt(sum(by15.whx(2:end).^2))), '  m']))
disp(strcat(['Water Temperature: ', num2str(sqrt(sum(by15.wtempx(2:end).^2))), '  C']))
disp(strcat(['Air Temperature: ', num2str(sqrt(sum(by15.atempx(2:end).^2))), '  C']))

disp('______________')
disp('2016')
disp('______________')
disp('--Annual Fit Means--')
disp('Wind Speed: 6.39 m/s')
disp('Wave Height: 2.27 m')
disp('Water Temperature: 16.84 C')
disp('Air Temperature: 15.80 C')
disp('--Annual Cycle Amplitude--')
disp(strcat(['Wind Speed: ', num2str(sqrt(sum(by16.wspdx(2:end).^2))), '  m/s']))
disp(strcat(['Wave Height: ', num2str(sqrt(sum(by16.whx(2:end).^2))), '  m']))
disp(strcat(['Water Temperature: ', num2str(sqrt(sum(by16.wtempx(2:end).^2))), '  C']))
disp(strcat(['Air Temperature: ', num2str(sqrt(sum(by16.atempx(2:end).^2))), '  C']))

%% Least Squares Fit with Semi-annual cycle

%least squares fit the four datasets for each year with a mean and annual cycle
f=[1/365,2/365];
%2015
[by15.wspdA2,by15.wspdx2]=LSfit(by15.dnum,by15.wspd,f);
[by15.whA2,by15.whx2]=LSfit(by15.dnum,by15.wh,f);
[by15.wtempA2,by15.wtempx2]=LSfit(by15.dnum,by15.wtemp,f);
[by15.atempA2,by15.atempx2]=LSfit(by15.dnum,by15.atemp,f);
%2016
[by16.wspdA2,by16.wspdx2]=LSfit(by16.dnum,by16.wspd,f);
[by16.whA2,by16.whx2]=LSfit(by16.dnum,by16.wh,f);
[by16.wtempA2,by16.wtempx2]=LSfit(by16.dnum,by16.wtemp,f);
[by16.atempA2,by16.atempx2]=LSfit(by16.dnum,by16.atemp,f);

%plot data along with LSfits
figure()
subplot(4,2,1)
plot(by15.dnum,by15.wspd,by15.dnum,by15.wspdA2*by15.wspdx2)
datetick('x','mmm')
title('2015 Semi-Annual Fitted Average Wind Speed')
ylabel('[m/s]')
subplot(4,2,3)
plot(by15.dnum,by15.wh,by15.dnum,by15.whA2*by15.whx2)
datetick('x','mmm')
title('"" Wave Height')
ylabel('[m]')
subplot(4,2,5)
plot(by15.dnum,by15.wtemp,by15.dnum,by15.wtempA2*by15.wtempx2)
datetick('x','mmm')
title('"" Water Temperature')
ylabel('[degC]')
subplot(4,2,7)
plot(by15.dnum,by15.atemp,by15.dnum,by15.atempA2*by15.atempx2)
datetick('x','mmm')
title(' "" Air Temperature')
ylabel('[degC]')


subplot(4,2,2)
plot(by16.dnum,by16.wspd,by16.dnum,by16.wspdA2*by16.wspdx2)
datetick('x','mmm')
title('2016 Semi-Annual Fitted Average Wind Speed')
ylabel('[m/s]')
subplot(4,2,4)
plot(by16.dnum,by16.wh,by16.dnum,by16.whA2*by16.whx2)
datetick('x','mmm')
title('"" Wave Height')
ylabel('[m]')
subplot(4,2,6)
plot(by16.dnum,by16.wtemp,by16.dnum,by16.wtempA2*by16.wtempx2)
datetick('x','mmm')
title(' "" Water Temperature')
ylabel('[degC]')
subplot(4,2,8)
plot(by16.dnum,by16.atemp,by16.dnum,by16.atempA2*by16.atempx2)
datetick('x','mmm')
title(' "" Air Temperature')
ylabel('[degC]')

% display mean and amplitude of the annual cycle for each year
disp('______________')
disp('2015')
disp('______________')
disp('--Semi-Annual Fit Means--')
disp('Wind Speed: 6.14 m/s')
disp('Wave Height: 2.00 m')
disp('Water Temperature: 17.94 C')
disp('Air Temperature: 16.79 C')
disp('--Annual Cycle Amplitude--')
disp(strcat(['Wind Speed: ', num2str(sqrt(sum(by15.wspdx2(2:3).^2))), 'm/s']))
disp(strcat(['Wave Height: ', num2str(sqrt(sum(by15.whx2(2:3).^2))), 'm']))
disp(strcat(['Water Temperature: ', num2str(sqrt(sum(by15.wtempx2(2:3).^2))), 'C']))
disp(strcat(['Air Temperature: ', num2str(sqrt(sum(by15.atempx2(2:3).^2))), 'C']))
disp('--Semi-Annual Cycle Amplitude--')
disp(strcat(['Wind Speed: ', num2str(sqrt(sum(by15.wspdx2(4:end).^2))), 'm/s']))
disp(strcat(['Wave Height: ', num2str(sqrt(sum(by15.whx2(4:end).^2))), 'm']))
disp(strcat(['Water Temperature: ', num2str(sqrt(sum(by15.wtempx2(4:end).^2))), 'C']))
disp(strcat(['Air Temperature: ', num2str(sqrt(sum(by15.atempx2(4:end).^2))), 'C']))

disp('______________')
disp('2016')
disp('______________')
disp('--Semi-Annual Fit Means--')
disp('Wind Speed: 6.36 m/s')
disp('Wave Height: 2.26 m')
disp('Water Temperature: 16.83 C')
disp('Air Temperature: 15.81 C')
disp('--Annual Cycle Amplitude--')
disp(strcat(['Wind Speed: ', num2str(sqrt(sum(by16.wspdx2(2:3).^2))), 'm/s']))
disp(strcat(['Wave Height: ', num2str(sqrt(sum(by16.whx2(2:3).^2))), 'm']))
disp(strcat(['Water Temperature: ', num2str(sqrt(sum(by16.wtempx2(2:3).^2))), 'C']))
disp(strcat(['Air Temperature: ', num2str(sqrt(sum(by16.atempx2(2:3).^2))), 'C']))
disp('--Semi-Annual Cycle Amplitude--')
disp(strcat(['Wind Speed: ', num2str(sqrt(sum(by16.wspdx2(4:end).^2))), 'm/s']))
disp(strcat(['Wave Height: ', num2str(sqrt(sum(by16.whx2(4:end).^2))), 'm']))
disp(strcat(['Water Temperature: ', num2str(sqrt(sum(by16.wtempx2(4:end).^2))), 'C']))
disp(strcat(['Air Temperature: ', num2str(sqrt(sum(by16.atempx2(4:end).^2))), 'C']))

%% Chi^2 and the Misfit


% annual cycle - 2015
by15.wspdchi2 = ((by15.wspd - (by15.wspdA*by15.wspdx)').^2) ./by15.wspdErr.^2;
by15.whchi2 = ((by15.wh - (by15.whA*by15.whx)').^2) ./by15.whErr.^2;
by15.wtempchi2 = ((by15.wtemp - (by15.wtempA*by15.wtempx)').^2) ./by15.wtempErr.^2;
by15.atempchi2 = ((by15.atemp - (by15.atempA*by15.atempx)').^2) ./by15.atempErr.^2;

% annual cycle - 2016
by16.wspdchi2 = ((by16.wspd - (by16.wspdA*by16.wspdx)').^2) ./by16.wspdErr.^2;
by16.whchi2 = ((by16.wh - (by16.whA*by16.whx)').^2) ./by16.whErr.^2;
by16.wtempchi2 = ((by16.wtemp - (by16.wtempA*by16.wtempx)').^2) ./by16.wtempErr.^2;
by16.atempchi2 = ((by16.atemp - (by16.atempA*by16.atempx)').^2) ./by16.atempErr.^2;

% semi-annual cycle - 2015
by15.wspd2chi2 = ((by15.wspd - (by15.wspdA2*by15.wspdx2)').^2) ./by15.wspdErr.^2;
by15.wh2chi2 = ((by15.wh - (by15.whA2*by15.whx2)').^2) ./by15.whErr.^2;
by15.wtemp2chi2 = ((by15.wtemp - (by15.wtempA2*by15.wtempx2)').^2) ./by15.wtempErr.^2;
by15.atemp2chi2 = ((by15.atemp - (by15.atempA2*by15.atempx2)').^2) ./by15.atempErr.^2;

% semi-annual cycle - 2016
by16.wspd2chi2 = ((by16.wspd - (by16.wspdA2*by16.wspdx2)').^2) ./by16.wspdErr.^2;
by16.wh2chi2 = ((by16.wh - (by16.whA2*by16.whx2)').^2) ./by16.whErr.^2;
by16.wtemp2chi2 = ((by16.wtemp - (by16.wtempA2*by16.wtempx2)').^2) ./by16.wtempErr.^2;
by16.atemp2chi2 = ((by16.atemp - (by16.atempA2*by16.atempx2)').^2) ./by16.atempErr.^2;

%plot chi^2 for annual and semi-annual cycles
figure()
subplot(4,1,1)
plot([2015, 2016],[sum(by15.wspdchi2), sum(by16.wspdchi2)], 'bo','LineWidth', 2)
hold on
plot([2015, 2016],[sum(by15.wspd2chi2), sum(by16.wspd2chi2)], 'rp','LineWidth', 2)
legend('Annual Fit','Semi-Annual Fit')
title('Wind Speed Misfit')
xlims([2014 2017])
ylabel('Misfit (m/s)')
xlabel('Year')

subplot(4,1,2)
plot([2015, 2016],[sum(by15.whchi2), sum(by16.whchi2)], 'bo','LineWidth', 2)
hold on
plot([2015, 2016],[sum(by15.wh2chi2), sum(by16.wh2chi2)], 'rp','LineWidth', 2)
title('Wave Height Misfit')
xlims([2014 2017])
ylabel('Misfit (m)')
xlabel('Year')

subplot(4,1,3)
plot([2015, 2016],[sum(by15.wtempchi2), sum(by16.wtempchi2)], 'bo','LineWidth', 2)
hold on
plot([2015, 2016],[sum(by15.wtemp2chi2), sum(by16.wtemp2chi2)], 'rp','LineWidth', 2)
title('Water Temperature Misfit')
xlims([2014 2017])
ylabel('Misfit (degC)')
xlabel('Year')

subplot(4,1,4)
plot([2015, 2016],[sum(by15.atempchi2), sum(by16.atempchi2)], 'bo','LineWidth', 2)
hold on
plot([2015, 2016],[sum(by15.atemp2chi2), sum(by16.atemp2chi2)], 'rp','LineWidth', 2)
title('Air Temperature Misfit')
xlims([2014 2017])
ylabel('Misfit (degC)')
xlabel('Year')

% plot chi^2 distributions
figure()
subplot(4,2,1)
plot(by15.dnum,by15.wspdchi2,'b',by15.dnum,by15.wspd2chi2,'r')
legend('Annual Fit', 'Semiannual Fit')
datetick('x','mmm')
title('2015 Wind Speed Misfit')
ylabel('[m/s]')
subplot(4,2,3)
plot(by15.dnum,by15.whchi2,'b',by15.dnum,by15.wh2chi2,'r')
datetick('x','mmm')
title('2015 Wave Height Misfit')
ylabel('[m]')
subplot(4,2,5)
plot(by15.dnum,by15.wtempchi2,'b',by15.dnum,by15.wtemp2chi2,'r')
datetick('x','mmm')
title('2015 Water Temperature Misfit')
ylabel('[degC]')
subplot(4,2,7)
plot(by15.dnum,by15.atempchi2,'b',by15.dnum,by15.atemp2chi2,'r')
datetick('x','mmm')
title('2015 Air Temperature Misfit')
ylabel('[degC]')


subplot(4,2,2)
plot(by16.dnum,by16.wspdchi2,'b',by16.dnum,by16.wspd2chi2,'r')
datetick('x','mmm')
title('2016 Wind Speed Misfit')
ylabel('[m/s]')
subplot(4,2,4)
plot(by16.dnum,by16.whchi2,'b',by16.dnum,by16.wh2chi2,'r')
datetick('x','mmm')
title('2016 Wave Height Misfit')
ylabel('[m]')
subplot(4,2,6)
plot(by16.dnum,by16.wtempchi2,'b',by16.dnum,by16.wtemp2chi2,'r')
datetick('x','mmm')
title('2016 Water Temperature Misfit')
ylabel('[degC]')
subplot(4,2,8)
plot(by16.dnum,by16.atempchi2,'b',by16.dnum,by16.atemp2chi2,'r')
datetick('x','mmm')
title('2016 Air Temperature Misfit')
ylabel('[degC]')


