clear all
%% problem #1 -- statistics of sea surface temperature
%read in data for 2021
file='/Users/kylie/Downloads/Homework/221_Data/HW1/scripps_pier-2021.nc';
time_2021=ncread(file,'time'); %from meta: time given in seconds from 1/1/1970
    date0=datenum(1970,1,1);
    dnum_2021=double(time_2021)/3600/24 + date0;
temperature_2021=ncread(file,'temperature');
%temperature_2021(find(temperature_2021>50))=NaN;%filter out temp values >50degC

%plot timeseries
figure();clf; 
plot(dnum_2021, temperature_2021);
xlabel('Months of 2021');ylabel('Temperature (^oC)');
datetick('x','mmm');
set(gca,'FontSize',16);
%axis([min(dnum_2021) max(dnum_2021) 0 25]);

%calculate mean, standard deviation
mean_temp=mean(temperature_2021,'omitnan')%15.2680
std_temp=std(temperature_2021,'omitnan')%1.8277

%empirically compute probability density function
%const=1/(std_temp.*sqrt(2*pi));
%pdf=const*exp(-((temperature_2021-mean_temp).^2)/(2*(std_temp.^2)))+mean_temp;
[number_per_bins,bins]=hist(temperature_2021,1:40);
figure();clf;plot(bins,number_per_bins/sum(number_per_bins));
    xlabel('Temperature (^oC)');ylabel('Probability Density');
    set(gca,'FontSize',16);

%% problem 2 -- extending the record
%grab data from 2005-2021
years=[2005:2021];
flist={};%zeros(length(years),1);
for n=1:length(years)
    flist{n}=['/Users/kylie/Downloads/Homework/221_Data/HW1/scripps_pier-' num2str(years(n)) '.nc'];
end

time={};temperature={};
for n=1:length(years)
    time{n}=ncread(flist{n},'time');
        dnum{n}=double(time{n})/3600/24 + date0;
    temperature{n}=ncread(flist{n},'temperature');
        temperature{n}(temperature{n}>50)=NaN;%filter out temp values > 50degc
end

%concatenate all cells, there is probably a much better way to do this, but I don't know it
total_time=cat(1,dnum{1},dnum{2},dnum{3},dnum{4},dnum{5},dnum{6},dnum{7},dnum{8},dnum{9},dnum{10},dnum{11},dnum{12},dnum{13},dnum{14},dnum{15},dnum{16},dnum{17});
total_temp=cat(1,temperature{1},temperature{2},temperature{3},temperature{4},temperature{5},temperature{6},temperature{7},temperature{8},temperature{9},temperature{10},temperature{11},temperature{12},temperature{13},temperature{14},temperature{15},temperature{16},temperature{17});
    
%plot timeseries
figure();clf; 
plot(total_time, total_temp);
xlabel('Time');ylabel('Temperature (^oC)');
datetick('x','yy');
set(gca,'FontSize',16);    

%calculate mean, standard deviation
mean_temp_total=mean(total_temp,'omitnan')
std_temp_total=std(total_temp,'omitnan')

%empirically compute probability density function
[number_per_bins,bins]=hist(total_temp,1:40);
figure();clf;plot(bins,number_per_bins/sum(number_per_bins));
    xlabel('Temperature (^oC)');ylabel('Probability Density');
    set(gca,'FontSize',16);
