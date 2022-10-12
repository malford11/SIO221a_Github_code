% 1. Statistics of sea surface temperature. Download the 2021 sea surface temperature data for the Scripps Pier from the SCCOOS web site: http://sccoos.org/thredds/catalog/autoss/catalog.html. 
% Read the temperature data, and produce a line plot of the 2021 temperatures with appropriately labeled axes. What do you observe in this plot? 
% Compute the mean and standard deviation for the sea surface temperature data. What do these statistics tell you about the temperature in 2021? 
% Compute an empirical probability density function for sea surface temperature. Does it look like any of the distributions that we discussed in class? 

t = ncread("autoss_scripps_pier-2021.nc",'time');
temp = ncread("autoss_scripps_pier-2021.nc",'temperature');
df = struct;
df.t = t;
df.temp = temp;
date0=datenum(1970,1,1);
dnum=double(df.t)/3600/24+date0;
df.t_dnum = dnum;
idx = find(df.t_dnum >= datenum(2021,1,1) & df.t_dnum <= datenum(2021,12,12));
vars = {'t','temperature'};
df.Temp2021 = table(df.t_dnum(idx), df.temp(idx));
df.Temp2021.Properties.VariableNames = vars;

figure(1)
clf
plot(df.Temp2021.t,df.Temp2021.temperature)
datetick
ylabel('temperature (^oC)')
title('Scripps Pier Temperature 2021')

% Temperature rises steadily from January to August, with a rapid flucatuation in late April, and sharply declines in late july/early august. 
% The temperature rises again in August before decreasing again until January.

mean_2021 = nanmean(df.Temp2021.temperature)
std_2021 = nanstd(df.Temp2021.temperature)

% The mean temperature in 2021 was 17.3883 C and a standard deviation is +/- 2.4467 C of this value. 
% The mean and standard deviation does not tell us everything in this data set because there is not a normal distribution. 
% The mean value is not actually the most common temperature, it is only the average of all temperatures. 

%PDF
clf
figure(2)
tiledlayout(1,2)
nexttile
histogram(df.Temp2021.temperature,20)
title('Histogram')
nexttile
histogram(df.Temp2021.temperature,20,"Normalization","pdf");
title('PDF')

% This resembles a bimodal distribution
% 
% 2. Extending the record. Now extend your record for the temperature from 2005- 2021 and repeat the calculations from the first exercise. 
% (This is a good time to practice using a loop to go through each of the data files.) What do you observe in these results? 
% In what ways are the 2021 results different from the 2005-2021 results? Is 2021 unusual? Is the sharp temperature change in August 2016 unusual? 
% I observed an extreme value recored early in the record which I decided to remove from the final plot as it is unlikely to be a real value.
% These results show a similar annual cycle to that observed in 2021
% 2021 did not have temperatures as warm as several years prior. It also did not reach minimum temperatures as low as years prior as well.
% The mean temperature in 2021 was about .5 C lower than the entire time series' average and the standard deviation was approximately .3 C lower than the overall standard deviation.
% 2021 had more temperature values in the 14.5 C range than would have been expected by a normal distribution. 
% 2021 was not too unusual of a year because its mean temperature was within one standard deviation of the full time series mean.
% The sharp change appears normal in August 2016. I added a running mean to the line plot which helps to show the drop in temperature was on par with other years.

%Complete Scripps Pier temperature record from https://data.caloos.org/#metadata/110895/station/data
df.t = ncread("scripps-pier-automated-shore-sta_1dd8_b354_6796.nc",'time');
df.temp = ncread("scripps-pier-automated-shore-sta_1dd8_b354_6796.nc", 'sea_water_temperature');
df.temp_qc = ncread("scripps-pier-automated-shore-sta_1dd8_b354_6796.nc", 'sea_water_temperature_qc_agg');
dnum=double(df.t)/3600/24+date0;
df.t_dnum = dnum;
idx = find(df.t_dnum >= datenum(2005,1,1) & df.t_dnum <= datenum(2021,12,12));
vars = {'t','temperature', 'temperature_qc'};
df.TempFull = table(df.t_dnum(idx), df.temp(idx), df.temp_qc(idx));
df.TempFull.Properties.VariableNames = vars;

ncdisp("scripps-pier-automated-shore-sta_f6dc_49c5_7e47.nc")

figure(3)
clf
plot(df.TempFull.t,df.TempFull.temperature)
datetick
ylabel('temperature (^oC)')
title('Scripps Pier Temperature 2005-2021')

%Maybe want to do some QC
mean_Full = nanmean(df.TempFull.temperature)
std_Full = nanstd(df.TempFull.temperature)
idx = find(df.TempFull.temperature >= mean_Full + 3.5*std_Full); %Find extreme outliers
df.TempFull(idx,:) = [];

%Recalculate mean and std now without extreme outliers
mean_Full = nanmean(df.TempFull.temperature)
std_Full = nanstd(df.TempFull.temperature)

figure(4)
clf
plot(df.TempFull.t,df.TempFull.temperature)
hold on
rmean = movmean(df.TempFull.temperature,(60*60*12));%seconds in half a day
plot(df.TempFull.t,rmean,'color','black',"LineWidth",2)
xticks(datenum(2005:2022,1,1))
xticklabels(2005:2022)
xtickangle(45)
xlim([datenum(2005,1,1) datenum(2022,1,1)])
ylabel('temperature (^oC)')
title('Scripps Pier Temperature 2005-2021')
legend('temperature','12 hour running mean')

%PDF
clf
figure(5)
tiledlayout(1,2)
nexttile
histogram(df.TempFull.temperature,20)
title('Histogram')
nexttile
histogram(df.TempFull.temperature,20,'Normalization',"pdf")
title('PDF')

% This is a normal distribution with a slight left skew

%Is it normal in August?
mean(df.TempFull.temperature(datenum(2016,8,1):datenum(2016,9,1)))
std(df.TempFull.temperature(datenum(2016,8,1):datenum(2016,9,1)))
df.TempFull.datetime =  datevec(df.TempFull.t);
idx = find(df.TempFull.datetime(:,2) == 8 );
mean(df.TempFull.temperature(idx))
std(df.TempFull.temperature(idx))
figure(6)
clf
scatter(df.TempFull.t(idx),df.TempFull.temperature(idx))
datetick('x',10)
