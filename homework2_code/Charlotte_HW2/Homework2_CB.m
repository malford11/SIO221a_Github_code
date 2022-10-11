%%%%SIOC 221A Data Analysis
%%%%Homework 2
%%%%Charlotte Bellerjeau

%% Choose any two years to compare from the SCOOS pier data
pier10.temp=ncread('PierData\scripps_pier-2010.nc','temperature');
pier10.pres=ncread('PierData\scripps_pier-2010.nc','pressure');
pier10.time=ncread('PierData\scripps_pier-2010.nc','time');
pier10.dnum=double(pier10.time)/3600/24+datenum(1970,1,1);

pier20.temp=ncread('PierData\scripps_pier-2020.nc','temperature');
pier20.pres=ncread('PierData\scripps_pier-2020.nc','pressure');
pier20.time=ncread('PierData\scripps_pier-2020.nc','time');
pier20.dnum=double(pier20.time)/3600/24+datenum(1970,1,1);

% remove outliers
pier20.temp(pier20.temp>35) = pier20.temp(find(pier20.temp>35)-1);
pier10.temp(pier10.temp>35) = pier10.temp(find(pier10.temp>35)-1);

%% Plot temperature and pressure for both years
figure()
subplot(2,2,1)
plot(pier10.dnum,pier10.temp)
datetick('x','mmm')
ylabel('Temperature [C]')
title('2010 Temperature')
subplot(2,2,2)
plot(pier20.dnum,pier20.temp)
datetick('x','mmm')
ylabel('Temperature [C]')
title('2020 Temperature')
subplot(2,2,3)
plot(pier10.dnum,pier10.pres)
datetick('x','mmm')
ylabel('Pressure [dbar]')
title('2010 Pressure')
subplot(2,2,4)
plot(pier20.dnum,pier20.pres)
datetick('x','mmm')
ylabel('Pressure [dbar]')
title('2020 Pressure')


%% Calculate mean, variance, and standard error of the mean and variance for both years
Tstats20 = stats_CB(pier20.temp);
pstats20 = stats_CB(pier20.pres);
Tstats10 = stats_CB(pier10.temp);
pstats10 = stats_CB(pier10.pres);

disp(strcat('2010 mean temperature: ',num2str(Tstats10(1)),'C +/- ',num2str(Tstats10(2)),'C'))
disp(strcat('2020 mean temperature: ',num2str(Tstats20(1)),'C +/- ',num2str(Tstats20(2)),'C'))
disp(strcat('2010 mean pressure: ',num2str(pstats10(1)),'dbar +/- ',num2str(pstats10(2)),'dbar'))
disp(strcat('2020 mean pressure: ',num2str(pstats20(1)),'dbar +/- ',num2str(pstats20(2)),'dbar'))

%is the mean consistant between errorbars
figure()
yyaxis left
errorbar([2010, 2020], [Tstats10(1), Tstats20(1)],  [Tstats10(2), Tstats20(2)],'.')
ylabel('Mean Temperature')
yyaxis right
errorbar([2010, 2020], [pstats10(1), pstats20(1)],  [pstats10(2), pstats20(2)],'.')
ylabel('Mean Pressure')
xlims([2005 2025])
title('Means')
xlabel('Year')


%is the variance consistant between errorbars
figure()
yyaxis left
errorbar([2010, 2020], [Tstats10(3), Tstats20(3)],  [Tstats10(4), Tstats20(4)],'.')
ylabel('Temperature Variance')
yyaxis right
xlims([2005 2025])
errorbar([2010, 2020], [pstats10(3), pstats20(3)],  [pstats10(4), pstats20(4)],'.')
ylabel('Pressure Variance')
title('Variance')
xlabel('Year')


%% Subsample the data at once per day and repeat the above
sr10 = round(length(pier10.temp)/365); % sample rate - samples per day
subind10 = sr10.*(1:364);%indices at which to subsample (once per day)
pier10.subtemp = pier10.temp(subind10); %subsample temperature
pier10.subpres = pier10.pres(subind10); %subsample pressure

sr20 = round(length(pier20.temp)/365); % sample rate - samples per day
subind20 = sr20.*(1:364); %indices at which to subsample (once per day)
pier20.subtemp = pier20.temp(subind20); %subsample temperature
pier20.subpres = pier20.pres(subind20); %subsample pressure

Tstats20sub = stats_CB(pier20.subtemp);
pstats20sub = stats_CB(pier20.subpres);
Tstats10sub = stats_CB(pier10.subtemp);
pstats10sub = stats_CB(pier10.subpres);

disp('When subsampling at once per day...')
disp(strcat('2010 mean temperature: ',num2str(Tstats10sub(1)),'C +/- ',num2str(Tstats10sub(2)),'C'))
disp(strcat('2020 mean temperature: ',num2str(Tstats20sub(1)),'C +/- ',num2str(Tstats20sub(2)),'C'))
disp(strcat('2010 mean pressure: ',num2str(pstats10sub(1)),'dbar +/- ',num2str(pstats10sub(2)),'dbar'))
disp(strcat('2020 mean pressure: ',num2str(pstats20sub(1)),'dbar +/- ',num2str(pstats20sub(2)),'dbar'))

%is the subsampled mean consistant between errorbars
figure()
yyaxis left
errorbar([2010, 2020], [Tstats10sub(1), Tstats20sub(1)],  [Tstats10sub(2), Tstats20sub(2)],'.')
ylabel('Mean Temperature')
yyaxis right
errorbar([2010, 2020], [pstats10sub(1), pstats20sub(1)],  [pstats10sub(2), pstats20sub(2)],'.')
ylabel('Mean Pressure')
xlims([2005 2025])
title('Means')
xlabel('Year')

%is the subsampled variance consistant between errorbars
figure()
yyaxis left
errorbar([2010, 2020], [Tstats10sub(3), Tstats20sub(3)],  [Tstats10sub(4), Tstats20sub(4)],'.')
ylabel('Temperature Variance')
yyaxis right
xlims([2005 2025])
errorbar([2010, 2020], [pstats10sub(3), pstats20sub(3)],  [pstats10sub(4), pstats20sub(4)],'.')
ylabel('Pressure Variance')
title('Variance')
xlabel('Year')


%% Extreme Values
%create a PDFand CDF for each year with a lot of bins (n=1000)
[numbinT10,binsT10]=hist(pier10.temp,1000);
[numbinT20,binsT20]=hist(pier20.temp,1000);
PDFT10 = numbinT10./sum(numbinT10);
PDFT20 = numbinT20./sum(numbinT20);
CDFT10 = cumtrapz(PDFT10);
CDFT20 = cumtrapz(PDFT20);

%observed magnitude of a 2sigma temperature event
sig210 = Tstats10(1)+(2*sqrt(Tstats10(3)));
sig220 = Tstats20(1)+(2*sqrt(Tstats20(3)));

%observed probability of a 2sigma temperature event
Oprobsig210 = 1-CDFT10(find(binsT10>sig210,1,'first'));
Oprobsig220 = 1-CDFT20(find(binsT20>sig220,1,'first'));

%display likelihood of 3 sigma event given observed distribution
disp('Using the observed temperature distribution, the likelihood of an extreme temperature event 2 standard deviations outside the mean is:')
disp(strcat('From 2010 Data: ',num2str(Oprobsig210)))
disp(strcat('From 2020 Data: ',num2str(Oprobsig220)))

%display likelihood of 3 sigma event given observed distribution
disp('Using a Gaussian temperature distribution, the likelihood of an extreme temperature event 2 standard deviations outside the mean is:')
disp('From any Data: 0.05')


%% Observed PDFs and CDFs
%observed PDF - temperature
[numbinT10,binsT10]=hist(pier10.temp,15);
[numbinT20,binsT20]=hist(pier20.temp,15);
PDFT10 = numbinT10./sum(numbinT10);
PDFT20 = numbinT20./sum(numbinT20);

%observed PDF - pressure
[numbinp10,binsp10]=hist(pier10.pres,10);
[numbinp20,binsp20]=hist(pier20.pres,10);
PDFp10 = numbinp10./sum(numbinp10);
PDFp20 = numbinp20./sum(numbinp20);


%% Gaussian PDFs and CDFs with same mean and std 
%gaussian PDF - temperature
PDFgT10 = 1/(sqrt(Tstats10(3))*sqrt(2*pi)).*exp(-(binsT10-Tstats10(1)).^2./(2*sqrt(Tstats10(3))^2));
PDFgT20 = 1/(sqrt(Tstats20(3))*sqrt(2*pi)).*exp(-(binsT20-Tstats20(1)).^2./(2*sqrt(Tstats20(3))^2));

%gaussian PDF - pressure
binsgp20 = 2:0.1:5;
binsgp10 = 3.5:0.1:6.5;
PDFgp10 = 1/(sqrt(pstats10(3)*2*pi)).*exp(-(binsgp10-pstats10(1)).^2./(2*pstats10(3)));
PDFgp20 = 1/(sqrt(pstats20(3)*2*pi)).*exp(-(binsgp20-pstats20(1)).^2./(2*pstats20(3)));

%% Uniform PDF with the same mean and standard deviation
%generate random data for uniform distributions with observed mean and std
unidT10 = (rand(100000,1)-.5)*sqrt(Tstats10(3))/std(rand(100000,1))+Tstats10(1);
unidT20 = (rand(100000,1)-.5)*sqrt(Tstats20(3))/std(rand(100000,1))+Tstats20(1);
unidp10 = (rand(100000,1)-.5)*sqrt(pstats10(3))/std(rand(100000,1))+pstats10(1);
unidp20 = (rand(100000,1)-.5)*sqrt(pstats20(3))/std(rand(100000,1))+pstats20(1);

%generate PDFs for uniform distribution
[numbinT10,binsuT10]=hist(unidT10,12:21);
[numbinT20,binsuT20]=hist(unidT20,12:23);
[numbinp10,binsup10]=hist(unidp10,3:0.5:7);
[numbinp20,binsup20]=hist(unidp20,2:0.5:5);
PDFuT10 = numbinT10./sum(numbinT10);
PDFuT20 = numbinT20./sum(numbinT20);
PDFup10 = numbinp10./sum(numbinp10);
PDFup20 = numbinp20./sum(numbinp20);


%% Plot the observed PDFs for both pressure and temperature with gaussian and uniform distributions overlayed
%plot Temperature distributions with gaussian and uniform distributions
figure()
plot(binsT10,PDFT10,'-b','Linewidth',2)
hold on
plot(binsT10,PDFgT10,'--b','LineWidth',1)
plot(binsuT10,PDFuT10,'--*b','LineWidth',1)
plot(binsT20,PDFT20,'-g','Linewidth',2)
plot(binsT20,PDFgT20,'--g','LineWidth',1)
plot(binsuT20,PDFuT20,'--*g','LineWidth',1)
title('Temperature PDFs')
xlabel('Temperature [C]')
ylabel('Probability')
legend('2010: Observed','2010: Gaussian','2010: Uniform','2020: Observed','2020: Gaussian','2020: Uniform')



%plot pressure distributions with gaussian and uniform distributions
figure()
plot(binsp10,PDFp10,'-b','Linewidth',2)
hold on
plot(binsgp10,PDFgp10,'--b','LineWidth',1)
plot(binsup10,PDFup10,'--*b','LineWidth',1)
plot(binsp20,PDFp20,'-g','Linewidth',2)
plot(binsgp20,PDFgp20,'--g','LineWidth',1)
plot(binsup20,PDFup20,'--*g','LineWidth',1)
title('Pressure PDFs')
xlabel('Pressure [dbar]')
ylabel('Probability')
legend('2010: Observed','2010: Gaussian','2010: Uniform','2020: Observed','2020: Gaussian','2020: Uniform')




%% Bonus...