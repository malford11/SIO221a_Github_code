%% Problem 1 
addpath(genpath('/Users/benjaminliu/Desktop/Oceanography/Course_material/Fall2022/data_analysis/homework'))

date0=datenum(1950,1,1);



info1= ncinfo('OS_T8S110W_DM134A-20150425_D_WIND_10min.nc');
%OS_T8S110W_DM183A-20160321_D_WIND_10min.nc
%OS_T8S110W_DM231A-20170606_D_WIND_10min.nc
data.yr15.time = ncread('OS_T8S110W_DM134A-20150425_D_WIND_10min.nc','TIME');
data.yr15.date = double(data.yr15.time)+date0; %convert time into date 
data.yr15.date = datetime(data.yr15.date,'ConvertFrom','datenum');
data.yr15.WSPD = ncread('OS_T8S110W_DM134A-20150425_D_WIND_10min.nc','WSPD');
data.yr15.UWND = ncread('OS_T8S110W_DM134A-20150425_D_WIND_10min.nc','UWND');
data.yr15.WSPD=data.yr15.WSPD';
data.yr15.UWND = data.yr15.UWND';
data.yr15.gap= diff(data.yr15.time);

data.yr16.time = ncread('OS_T8S110W_DM183A-20160321_D_WIND_10min.nc','TIME');
data.yr16.date = double(data.yr16.time)+date0; %convert time into date 
data.yr16.date = datetime(data.yr16.date,'ConvertFrom','datenum');
data.yr16.WSPD = ncread('OS_T8S110W_DM183A-20160321_D_WIND_10min.nc','WSPD');
data.yr16.UWND = ncread('OS_T8S110W_DM183A-20160321_D_WIND_10min.nc','UWND');
data.yr16.WSPD=data.yr16.WSPD';
data.yr16.UWND = data.yr16.UWND';
data.yr16.gap= diff(data.yr16.time);

data.yr17.time = ncread('OS_T8S110W_DM231A-20170606_D_WIND_10min.nc','TIME');
data.yr17.date = double(data.yr17.time)+date0; %convert time into date 
data.yr17.date = datetime(data.yr17.date,'ConvertFrom','datenum');
data.yr17.WSPD = ncread('OS_T8S110W_DM231A-20170606_D_WIND_10min.nc','WSPD');
data.yr17.UWND = ncread('OS_T8S110W_DM231A-20170606_D_WIND_10min.nc','UWND');
data.yr17.WSPD=data.yr17.WSPD';
data.yr17.UWND = data.yr17.UWND';
data.yr17.gap= diff(data.yr17.time);

data.sum.time =[data.yr15.time; data.yr16.time; data.yr17.time];
data.sum.date =[data.yr15.date; data.yr16.date; data.yr17.date];
data.sum.WSPD =[data.yr15.WSPD; data.yr16.WSPD; data.yr17.WSPD];
data.sum.UWND =[data.yr15.UWND; data.yr16.UWND; data.yr17.UWND];
data.sum.gap =diff(data.sum.time);




figure(1)
plot(data.sum.date,data.sum.WSPD)
xlabel('Date','FontSize',16)
ylabel('Total Wind Speed (m/s)','FontSize',16)
figure(2)
plot(data.sum.date,data.sum.UWND)
xlabel('Date','FontSize',16)
ylabel('Total Wind Speed (m/s)','FontSize',16)

figure(3)
plot(data.sum.gap)
xlabel('Date Points','FontSize',16)
ylabel('Time Interval (days)','FontSize',16)


nan_ind_15 = ~isnan(data.yr15.WSPD);
data.yr15.WSPD_interp = interp1(data.yr15.time(nan_ind_15),data.yr15.WSPD(nan_ind_15),data.yr15.time,'linear');
data.yr15.WSPD_interp1 = data.yr15.WSPD_interp(1:47581); %delete the last 24 nan in interpreted data in 2015
nan_ind_15_1 = ~isnan(data.yr15.UWND);
data.yr15.UWND_interp = interp1(data.yr15.time(nan_ind_15_1),data.yr15.WSPD(nan_ind_15_1),data.yr15.time,'linear');


nan_ind_16 = ~isnan(data.yr16.WSPD);
data.yr16.WSPD_interp = interp1(data.yr16.time(nan_ind_16),data.yr16.WSPD(nan_ind_16),data.yr16.time,'linear');
nan_ind_16_1 = ~isnan(data.yr16.UWND);
data.yr16.UWND_interp = interp1(data.yr16.time(nan_ind_16_1),data.yr16.WSPD(nan_ind_16_1),data.yr16.time,'linear');

nan_ind_17 = ~isnan(data.yr17.WSPD);
data.yr17.WSPD_interp = interp1(data.yr17.time(nan_ind_17),data.yr17.WSPD(nan_ind_17),data.yr17.time,'linear');
nan_ind_17_1 = ~isnan(data.yr17.UWND);
data.yr17.UWND_interp = interp1(data.yr17.time(nan_ind_17_1),data.yr17.WSPD(nan_ind_17_1),data.yr17.time,'linear');


data.sum.WSPD_interp =[data.yr15.WSPD_interp; data.yr16.WSPD_interp; data.yr17.WSPD_interp];
data.sum.UWND_interp =[data.yr15.UWND_interp; data.yr16.UWND_interp; data.yr17.UWND_interp];

%% Problem 2 



N = 8640/2; %8640:number of data points (sampling fre 10min) corresponding to 60 days)
%2015 data 
nbin_15 = floor((length(data.yr15.time)-N)/N);
for ii=1:nbin_15
    data.yr15.WSPD_window{ii}= data.yr15.WSPD_interp(((ii-1)*N+1):((ii-1)*N+1)+2*N);
    data.yr15.UWND_window{ii}= data.yr15.UWND_interp(((ii-1)*N+1):((ii-1)*N+1)+2*N);
end 
% data.yr15.WSPD_window{nbin_15+1}= data.yr15.WSPD((((ii-1)*N+1)+2*N)+1:length(data.yr15.time)); %take the rest of the data(less than 60days) as the last segment
% data.yr15.UWND_window{nbin_15+1}= data.yr15.UWND((((ii-1)*N+1)+2*N)+1:length(data.yr15.time));



nbin_16 = floor((length(data.yr16.time)-N)/N);
for ii=1:nbin_16
    data.yr16.WSPD_window{ii}= data.yr16.WSPD_interp(((ii-1)*N+1):((ii-1)*N+1)+2*N);
    data.yr16.UWND_window{ii}= data.yr16.UWND_interp(((ii-1)*N+1):((ii-1)*N+1)+2*N);
end 
% data.yr16.WSPD_window{nbin_16+1}= data.yr16.WSPD((((ii-1)*N+1)+2*N)+1:length(data.yr16.time)); %take the rest of the data(less than 60days) as the last segment
% data.yr16.UWND_window{nbin_16+1}= data.yr16.UWND((((ii-1)*N+1)+2*N)+1:length(data.yr16.time));



nbin_17 = floor((length(data.yr17.time)-N)/N);
for ii=1:nbin_17
    data.yr17.WSPD_window{ii}= data.yr17.WSPD_interp(((ii-1)*N+1):((ii-1)*N+1)+2*N);
    data.yr17.UWND_window{ii}= data.yr17.UWND_interp(((ii-1)*N+1):((ii-1)*N+1)+2*N);
end 
% data.yr17.WSPD_window{nbin_17+1}= data.yr17.WSPD((((ii-1)*N+1)+2*N)+1:length(data.yr17.time)); %take the rest of the data(less than 60days) as the last segment
% data.yr17.UWND_window{nbin_17+1}= data.yr17.UWND((((ii-1)*N+1)+2*N)+1:length(data.yr17.time));

% data.yr15.WSPD_window
% data.yr16.WSPD_window
% data.yr17.WSPD_window
% data.yr15.UWND_window
% data.yr16.UWND_window
% data.yr17.UWND_window
%% Problem 3 
dt = 0.006944444441615; % in days (10min)
for j = 1:length(data.yr15.WSPD_window)

    % raw data subtract mean 
data_demean{j} = data.yr15.WSPD_window{j}-mean(data.yr15.WSPD_window{j});

[f_demean(j,:),xs_demean(j,:),spec_demean(j,:)] = fft_data(data_demean{j},dt);
figure(4)
loglog(f_demean(j,:),spec_demean(j,:))
xlabel('Frequency (cycle per day)','FontSize',16)
ylabel('\Phi_v (m/s)^2/cpd','FontSize',16)
title('Raw Data Demean')

hold on 

%detrended data 
data_detre{j} = detrend(data.yr15.WSPD_window{j});
[f_detre(j,:),xs_detre(j,:),spec_detre(j,:)] = fft_data(data_detre{j},dt);

figure(5)
loglog(f_detre(j,:),spec_detre(j,:))
xlabel('Frequency (cycle per day)','FontSize',16)
ylabel('\Phi_v (m/s)^2/cpd','FontSize',16)
title('Detrended Data')

hold on 

%Hanning window applied 
%normf{j}=sqrt((sum(hann(length(data.yr15.WSPD_window{j})).^2))/length(data.yr15.WSPD_window{j}));
normf = sqrt(8/3);
data_hanning{j} = (hann(length(data.yr15.WSPD_window{j}))*normf).*detrend(data.yr15.WSPD_window{j});
[f_hanning(j,:),xs_hanning(j,:),spec_hanning(j,:)] = fft_data(data_hanning{j},dt);

figure(6)
loglog(f_hanning(j,:),spec_hanning(j,:))
xlabel('Frequency (cycle per day)','FontSize',16)
ylabel('\Phi_v (m/s)^2/cpd','FontSize',16)
title('Detrended Hanning Window Data')

hold on 

end 
spec_demean_ave = mean(spec_demean);
spec_detre_ave = mean(spec_detre);
spec_hanning_ave = mean(spec_hanning);

figure(7)
loglog(f_demean(1,:),spec_demean_ave)
xlabel('Frequency (cycle per day)','FontSize',16)
ylabel('\Phi_v (m/s)^2/cpd','FontSize',16)
title('Raw Data Demean')


figure(8)
loglog(f_detre(1,:),spec_detre_ave)
xlabel('Frequency (cycle per day)','FontSize',16)
ylabel('\Phi_v (m/s)^2/cpd','FontSize',16)
title('Detrended Data')


figure(9)
loglog(f_hanning(1,:),spec_hanning_ave)
xlabel('Frequency (cycle per day)','FontSize',16)
ylabel('\Phi_v (m/s)^2/cpd','FontSize',16)
title('Detrended Hanning Window Data')

%% Problem 4 

nu=2*nbin_15;
err_low = nu/chi2inv(1-.05/2,nu);
err_high = nu/chi2inv(.05/2,nu);
uncertainty_spec_hanning_ave(1,:) = err_high *spec_hanning_ave; % 1 raw is the lower bound 
uncertainty_spec_hanning_ave(2,:) = err_low *spec_hanning_ave; %2 raw is the upper bound 

figure(10)
errorbar(f_hanning(1,:),spec_hanning_ave,spec_hanning_ave-uncertainty_spec_hanning_ave(2,:),uncertainty_spec_hanning_ave(1,:)-spec_hanning_ave)
set(gca, 'XScale','log', 'YScale','log')
ylim([10^0,10^8])
xlabel('Frequency (cycle per day)','FontSize',16)
ylabel('\Phi_v (m/s)^2/cpd','FontSize',16)
title('Hanning Window Data')


uncertainty_spec_demean_ave(1,:) = err_high *spec_demean_ave; % 1 raw is the lower bound 
uncertainty_spec_demean_ave(2,:) = err_low *spec_demean_ave; %2 raw is the upper bound 

figure(11)
errorbar(f_demean(1,:),spec_demean_ave,spec_demean_ave-uncertainty_spec_demean_ave(2,:),uncertainty_spec_demean_ave(1,:)-spec_demean_ave)
set(gca, 'XScale','log', 'YScale','log')
ylim([10^0,10^8])
xlabel('Frequency (cycle per day)','FontSize',16)
ylabel('\Phi_v (m/s)^2/cpd','FontSize',16)
title('Demean Data')



uncertainty_spec_detre_ave(1,:) = err_high *spec_detre_ave; % 1 raw is the lower bound 
uncertainty_spec_detre_ave(2,:) = err_low *spec_detre_ave; %2 raw is the upper bound 

figure(12)
errorbar(f_detre(1,:),spec_detre_ave,spec_detre_ave-uncertainty_spec_detre_ave(2,:),uncertainty_spec_detre_ave(1,:)-spec_detre_ave)
set(gca, 'XScale','log', 'YScale','log')
ylim([10^0,10^8])
xlabel('Frequency (cycle per day)','FontSize',16)
ylabel('\Phi_v (m/s)^2/cpd','FontSize',16)
title('Detrended Data')


%% Problem 5

[f_16,spec_ave_16,uncertainty_16] = ave_fft_uncertainty(data.yr16.WSPD_window,dt,95);
[f_17,spec_ave_17,uncertainty_17] = ave_fft_uncertainty(data.yr17.WSPD_window,dt,95);

figure(13)
errorbar(f_hanning(1,:),spec_hanning_ave,spec_hanning_ave-uncertainty_spec_hanning_ave(2,:),uncertainty_spec_hanning_ave(1,:)-spec_hanning_ave)
hold on 
errorbar(f_16.hanning,spec_ave_16.hanning,spec_ave_16.hanning-uncertainty_16.hanning(2,:),uncertainty_16.hanning(1,:)-spec_ave_16.hanning)
hold on 
errorbar(f_17.hanning,spec_ave_17.hanning,spec_ave_17.hanning-uncertainty_17.hanning(2,:),uncertainty_17.hanning(1,:)-spec_ave_17.hanning)

set(gca, 'XScale','log', 'YScale','log')
%ylim([10^0,10^8])
legend('2015 data','2016 data','2017 data')
title('Hanning Window Data')

xlabel('Frequency (cycle per day)','FontSize',16)
ylabel('\Phi_v (m/s)^2/cpd','FontSize',16)

figure(14)
errorbar(f_detre(1,:),spec_detre_ave,spec_detre_ave-uncertainty_spec_detre_ave(2,:),uncertainty_spec_detre_ave(1,:)-spec_detre_ave)
hold on 
errorbar(f_16.detre,spec_ave_16.detre,spec_ave_16.detre-uncertainty_16.detre(2,:),uncertainty_16.detre(1,:)-spec_ave_16.detre)
hold on 
errorbar(f_17.detre,spec_ave_17.detre,spec_ave_17.detre-uncertainty_17.detre(2,:),uncertainty_17.detre(1,:)-spec_ave_17.detre)

set(gca, 'XScale','log', 'YScale','log')
ylim([10^0,10^6])
legend('2015 data','2016 data','2017 data')
title('Detrended Data')
xlabel('Frequency (cycle per day)','FontSize',16)
ylabel('\Phi_v (m/s)^2/cpd','FontSize',16)


% figure(15)
% plot(f_hanning(1,:),spec_hanning_ave)
% hold on 
% plot(f_16.hanning,spec_ave_16.hanning)
% hold on 
% plot(f_17.hanning,spec_ave_17.hanning)
% 
% set(gca, 'XScale','log', 'YScale','log')
% %ylim([10^0,10^8])
% legend('2015 data','2016 data','2017 data')
% xlabel('Frequency (cycle per day)','FontSize',16)
% ylabel('\Phi_v (m/s)^2/cpd','FontSize',16)


%% Problem 6
[f_15_UWND,spec_UWND_ave_15,uncertainty_UWND_15] = ave_fft_uncertainty(data.yr15.UWND_window,dt,95);
[f_16_UWND,spec_UWND_ave_16,uncertainty_UWND_16] = ave_fft_uncertainty(data.yr16.UWND_window,dt,95);
[f_17_UWND,spec_UWND_ave_17,uncertainty_UWND_17] = ave_fft_uncertainty(data.yr17.UWND_window,dt,95);


figure(16)
errorbar(f_15_UWND.detre,spec_UWND_ave_15.detre,spec_UWND_ave_15.detre-uncertainty_UWND_15.detre(2,:),uncertainty_UWND_15.detre(1,:)-spec_UWND_ave_15.detre)
hold on 
errorbar(f_16_UWND.detre,spec_UWND_ave_16.detre,spec_UWND_ave_16.detre-uncertainty_UWND_16.detre(2,:),uncertainty_UWND_16.detre(1,:)-spec_UWND_ave_16.detre)
hold on 
errorbar(f_17_UWND.detre,spec_UWND_ave_17.detre,spec_UWND_ave_17.detre-uncertainty_UWND_17.detre(2,:),uncertainty_UWND_17.detre(1,:)-spec_UWND_ave_17.detre)

set(gca, 'XScale','log', 'YScale','log')
ylim([10^0,10^6])
title('Detrended Data')
legend('2015 data','2016 data','2017 data')
xlabel('Frequency (cycle per day)','FontSize',16)
ylabel('\Phi_v (m/s)^2/cpd','FontSize',16)

figure(17)
errorbar(f_15_UWND.hanning,spec_UWND_ave_15.hanning,spec_UWND_ave_15.hanning-uncertainty_UWND_15.hanning(2,:),uncertainty_UWND_15.hanning(1,:)-spec_UWND_ave_15.hanning)
hold on 
errorbar(f_16_UWND.hanning,spec_UWND_ave_16.hanning,spec_UWND_ave_16.hanning-uncertainty_UWND_16.hanning(2,:),uncertainty_UWND_16.hanning(1,:)-spec_UWND_ave_16.hanning)
hold on 
errorbar(f_17_UWND.hanning,spec_UWND_ave_17.hanning,spec_UWND_ave_17.hanning-uncertainty_UWND_17.hanning(2,:),uncertainty_UWND_17.hanning(1,:)-spec_UWND_ave_17.hanning)

set(gca, 'XScale','log', 'YScale','log')
ylim([10^0,10^8])
title('Hanning Window Data')
legend('2015 data','2016 data','2017 data')
xlabel('Frequency (cycle per day)','FontSize',16)
ylabel('\Phi_v (m/s)^2/cpd','FontSize',16)


figure(18)
errorbar(f_15_UWND.detre,spec_UWND_ave_15.detre,spec_UWND_ave_15.detre-uncertainty_UWND_15.detre(2,:),uncertainty_UWND_15.detre(1,:)-spec_UWND_ave_15.detre,'LineWidth',3)
hold on 
errorbar(f_16_UWND.detre,spec_UWND_ave_16.detre,spec_UWND_ave_16.detre-uncertainty_UWND_16.detre(2,:),uncertainty_UWND_16.detre(1,:)-spec_UWND_ave_16.detre,'LineWidth',3)
hold on 
errorbar(f_17_UWND.detre,spec_UWND_ave_17.detre,spec_UWND_ave_17.detre-uncertainty_UWND_17.detre(2,:),uncertainty_UWND_17.detre(1,:)-spec_UWND_ave_17.detre,'LineWidth',3)

set(gca, 'YScale','log')
ylim([10^0,10^6])
title('Detrended Data')
legend('2015 data','2016 data','2017 data')
xlabel('Frequency (cycle per day)','FontSize',16)
ylabel('\Phi_v (m/s)^2/cpd','FontSize',16)
xlim([0.95,1.05])


figure(19)
errorbar(f_15_UWND.hanning,spec_UWND_ave_15.hanning,spec_UWND_ave_15.hanning-uncertainty_UWND_15.hanning(2,:),uncertainty_UWND_15.hanning(1,:)-spec_UWND_ave_15.hanning,'LineWidth',3)
hold on 
errorbar(f_16_UWND.hanning,spec_UWND_ave_16.hanning,spec_UWND_ave_16.hanning-uncertainty_UWND_16.hanning(2,:),uncertainty_UWND_16.hanning(1,:)-spec_UWND_ave_16.hanning,'LineWidth',3)
hold on 
errorbar(f_17_UWND.hanning,spec_UWND_ave_17.hanning,spec_UWND_ave_17.hanning-uncertainty_UWND_17.hanning(2,:),uncertainty_UWND_17.hanning(1,:)-spec_UWND_ave_17.hanning,'LineWidth',3)

set(gca, 'YScale','log')
ylim([10^0,10^8])
title('Hanning Window Data')
legend('2015 data','2016 data','2017 data')
xlabel('Frequency (cycle per day)','FontSize',16)
ylabel('\Phi_v (m/s)^2/cpd','FontSize',16)
xlim([0.95,1.05])

[test1,test2] = fft_ave_uncertainty(data.yr15.time,dt,data.yr15.WSPD_interp,N,95,3);