% script for final project for SIOC221A
% for each of four moorings plot
%   - frequency spectra with depth
%   - M2 peak with depth

%% Data Visualization

lon=ncread('blt_canyon_mb_qc.nc','lon');
lat=ncread('blt_canyon_mb_qc.nc','lat');
elv=-ncread('blt_canyon_mb_qc.nc','z');

moorlat = [54.197,54.204,54.185,54.190,54.239];
moorlon = [-11.862,-11.873,-11.846,-11.852,-11.949];

moorname = {'MAVS1','MP2','MP3','TCHAIN','MP1'};
moordates = {'7/6/21-10/9/21','7/7/21-10/6/21','10/21/21-8/2/22','7/6/21-8/11/22','6/28/21-7/5/21'};

figure()
contour(lon,lat,elv',-1000:-25:-3000,'-k','LineWidth',1,'Color',[0.5 0.5 0.5])
hold on
plot(moorlon(3),moorlat(3),'gd','LineWidth',2);
plot(moorlon(4),moorlat(4),'rp','LineWidth',2);
plot(moorlon(1:2),moorlat(1:2),'bo','LineWidth',2);
%plot(moorlon(5),moorlat(5),'yo','Linewidth',2);
text(moorlon(1)+0.004,moorlat(1),moorname(1),'Linewidth',2)
text(moorlon(2)+0.004,moorlat(2),moorname(2),'Linewidth',2)
text(moorlon(3)-0.013,moorlat(3)+0.0015,moorname(3),'Linewidth',2)
text(moorlon(4)+0.004,moorlat(4)+0.001,moorname(4),'Linewidth',2)
%text(moorlon(5)+0.004,moorlat(5),moorname(5),'Linewidth',2)
title('Data Visualization')
xlims([-11.98 -11.8])
ylims([54.155 54.25])
legend('-','10 Month','12 Month','3 Month')

% find depth at each mooring
for i=1:length(moorlon)
    xind = find(lon>moorlon(i),1,'first');
    yind = find(lat>moorlat(i),1,'first');
    moorbot(i)=abs(elv(xind,yind));
end

%% MP3
%load mooring data
trim = [1 10 27248; 1 2 35];
mp3 = mooring_load('blt_mp3.nc', trim);
mp3.time = mp3.time./1e6./60./60./24;  % time is now in days since 2021-10-15
mp3.z=ncread('blt_mp3.nc','depth');  %some mooring files have 'depth' while others have 'z'
mp3.z(1:trim(2,2))=[];
mp3.z(trim(2,3):end)=[];
% %add dnum time vector for mp3 mooring
% mp3.dnum = yday2datenum(mp3.time+295.3285,2021); % time vector in datenum format
[mmp3,nmp3] = size(mp3.u); %[time,depth]

%add WKB stretched and scaled components
mp3 = WKB_BLTstrat(mp3);

%add frequency spectra of upcanyon velocity at each depth
mp3 = mooring_spectra_221A(mp3,6000);
%spectral uncertainty
nu=2*8;
err_high = nu/chi2inv(.05/2,nu);
err_low =  nu/chi2inv(1-.05/2,nu);

% plot frequency spectra at each depth
DepthSpectra = cbrewer('div','RdYlGn',nmp3);
figure()
for i=1:nmp3
    plot3(mp3.k,mp3.z(i)*ones(length(mp3.k),1),mp3.Fs(:,i),'-','Color',DepthSpectra(i,:))
    hold on
end
plot3([2e-2 2e-2],[mean(mp3.z) mean(mp3.z)],[err_low err_high].*1.9,'k-*')
set(gca,'XScale','log')
set(gca,'ZScale','log')
title('MP3 Frequency Spectra')
xlabel('Frequency [cpd]')
ylabel('Depth [m]')
zlabel('\Phi_v (m/s)^2/cpd')

% find M2 peak over depth
mp3.M2peak = M2peak(mp3.Fs,mp3.k,nmp3);
mp3.zAB = moorbot(3)-mp3.z;

figure()
for j = 1:length(mp3.z)
    plot3(mp3.k,mp3.z(j)*ones(length(mp3.k),1),mp3.Fs(:,j),'Color',DepthSpectra(j,:))
    plot3(1.95,mp3.z(j),mp3.M2peak(j),'ko')
    hold on
end
set(gca,'XScale','log')
set(gca,'ZScale','log')
title('M2 Peaks on MP3 Spectra')
xlabel('Frequency [cpd]')
ylabel('Depth [m]')
zlabel('\Phi_v (m/s)^2/cpd')




%% TCHAIN
%load mooring data
trim = [1 48 38387; 1 3 7];
tchain = mooring_load('blt_tchain.nc', trim);
tchain.time = tchain.time./1e6./60./60./24;  % time is now in days since 2021-07-06 5am
tchain.z=ncread('blt_tchain.nc','depth');  %some mooring files have 'depth' while others have 'z'
tchain.z(1:trim(2,2))=[];
tchain.z(trim(2,3):end)=[];
% %add dnum time vector for mp3 mooring
% tchain.dnum = yday2datenum(tchain.time+187.2083,2021); % time vector in datenum format
[mtchain,ntchain] = size(tchain.u); %[time,depth]

%add WKB stretched and scaled components
tchain = WKB_BLTstrat(tchain);

%add frequency spectra of upcanyon velocity at each depth
tchain = mooring_spectra_221A(tchain,6000);
%spectral uncertainty
nu=2*12;
err_high = nu/chi2inv(.05/2,nu);
err_low =  nu/chi2inv(1-.05/2,nu);

% plot frequency spectra at each depth
DepthSpectra = cbrewer('div','RdYlGn',ntchain);
figure()
for i=1:ntchain
    %counterclockwise spectrum
    plot3(tchain.k,tchain.z(i)*ones(length(tchain.k),1),tchain.Fs(:,i),'-','Color',DepthSpectra(i,:))
    hold on
end
plot3([2e-2 2e-2],[mean(tchain.z) mean(tchain.z)],[err_low err_high].*1.9,'k-*')
set(gca,'XScale','log')
set(gca,'ZScale','log')
title('TCHAIN Frequency Spectra')
xlabel('Frequency [cpd]')
ylabel('Depth [m]')
zlabel('\Phi_v (m/s)^2/cpd')

% find M2 peak over depth
tchain.M2peak = M2peak(tchain.Fs,tchain.k,ntchain);
tchain.zAB = moorbot(4)-tchain.z;



%% MAVS1

%load mooring data
trim = [1 40 8717; 1 2 21];
mavs1 = mooring_load('blt_mavs1.nc', trim);
%time is in 'minutes since 2021-07-06 05:00:00'
mavs1.time = mavs1.time./60./24;  % time is now in days since 2021-07-06 05:00:00
mavs1.z=ncread('blt_mavs1.nc','z');  %some mooring files have 'depth' while others have 'z'
mavs1.z(1:trim(2,2))=[];
mavs1.z(trim(2,3):end)=[];
% %add dnum time vector for mp3 mooring
% mavs1.dnum = yday2datenum(mavs1.time+187.2083,2021); % time vector in datenum format
[mmavs1,nmavs1] = size(mavs1.u); %[time,depth]

%add WKB stretched and scaled components
mavs1 = WKB_BLTstrat(mavs1);

%add frequency and wavenumber spectra of velocity
mavs1 = mooring_spectra_221A(mavs1,2000);
%spectral uncertainty
nu=2*8;
err_high = nu/chi2inv(.05/2,nu);
err_low =  nu/chi2inv(1-.05/2,nu);


% plot frequency spectra at each depth
DepthSpectra = cbrewer('div','RdYlGn',nmavs1);
figure()
for i=1:nmavs1
    %counterclockwise spectrum
    plot3(mavs1.k,mavs1.z(i)*ones(length(mavs1.k),1),mavs1.Fs(:,i),'-','Color',DepthSpectra(i,:))
    hold on
end
plot3([2e-2 2e-2],[mean(mavs1.z) mean(mavs1.z)],[err_low err_high].*1.9,'k-*')
set(gca,'XScale','log')
set(gca,'ZScale','log')
title('MAVS1 Frequency Spectra')
xlabel('Frequency [cpd]')
ylabel('Depth [m]')
zlabel('\Phi_v (m/s)^2/cpd')

% find M2 peak over depth
mavs1.M2peak = M2peak(mavs1.Fs,mavs1.k,nmavs1);
mavs1.zAB = moorbot(1)-mavs1.z;



%% MP2

%load mooring data
trim = [1 27 8630; 1 6 34];
mp2 = mooring_load('blt_mp2.nc', trim);
[mmp2,nmp2] = size(mp2.u); %[time,depth]
%time is in 'minutes since 2021-07-07 10:15:00'
mp2.time = mp2.time./60./24;  % time is now in days since 2021-07-07 10:15:00
mp2.z=ncread('blt_mp2.nc','z');  %some mooring files have 'depth' while others have 'z'
mp2.z(1:trim(2,2))=[];
mp2.z(trim(2,3):end)=[];
% %add dnum time vector for mp3 mooring
% mp2.dnum = yday2datenum(mp2.time+188.4271,2021); % time vector in datenum format

%add WKB stretched and scaled components
mp2 = WKB_BLTstrat(mp2);

%add frequency and wavenumber spectra of velocity
mp2 = mooring_spectra_221A(mp2,2000);
%spectral uncertainty
nu=2*8;
err_high = nu/chi2inv(.05/2,nu);
err_low =  nu/chi2inv(1-.05/2,nu);

% plot frequency spectra at each depth
DepthSpectra = cbrewer('div','RdYlGn',nmp2);
figure()
for i=1:nmp2
    %counterclockwise spectrum
    plot3(mp2.k,mp2.z(i)*ones(length(mp2.k),1),mp2.Fs(:,i),'-','Color',DepthSpectra(i,:))
    hold on
end
plot3([2e-2 2e-2],[mean(mp2.z) mean(mp2.z)],[err_low err_high].*1.9,'k-*')
set(gca,'XScale','log')
set(gca,'ZScale','log')
title('MP2 Frequency Spectra')
xlabel('Frequency [cpd]')
ylabel('Depth [m]')
zlabel('\Phi_v (m/s)^2/cpd')

% find M2 peak over depth
mp2.M2peak = M2peak(mp2.Fs,mp2.k,nmp2);
mp2.zAB = moorbot(2)-mp2.z;

% %% MP1
% %load mooring data
% trim = [1 11 633; 1 4 38];
% mp1 = mooring_load('blt_mp1.nc', trim);
% [mmp1,nmp1] = size(mp1.u); %[time,depth]
% %time is in 'minutes since 2021-06-28 15:15:00
% mp1.time = mp1.time./60./24;  % time is now in days since 2021-06-28 15:15:00
% mp1.z=ncread('blt_mp1.nc','z');  %some mooring files have 'depth' while others have 'z'
% mp1.z(1:trim(2,2))=[];
% mp1.z(trim(2,3):end)=[];
% %add dnum time vector for mp3 mooring
% mp1.dnum = yday2datenum(mp1.time+179.6354,2021); % time vector in datenum format
% 
% %add WKB stretched and scaled components
% mp1 = WKB_BLTstrat(mp1);
% 
% %add frequency and wavenumber spectra of velocity
% mp1 = mooring_spectra_221A(mp1);
% 
% % plot frequency spectra at each depth
% DepthSpectra = cbrewer('div','RdYlGn',nmp1);
% figure()
% for i=1:nmp1
%     %counterclockwise spectrum
%     plot3(mp1.k,mp1.z(i)*ones(length(mp1.k),1),mp1.Fs(:,i),'-','Color',DepthSpectra(i,:))
%     hold on
% end
% set(gca,'XScale','log')
% set(gca,'ZScale','log')
% title('MP1 Frequency Spectra')
% xlabel('Frequency [cpd]')
% ylabel('Depth [m]')
% zlabel('\Phi_v (m/s)^2/cpd')
% 
% % find M2 peak over depth
% mp1.M2peak = M2peak(mp1.Fs,mp1.k,nmp1);
% mp1.zAB = moorbot(5)-mp1.z;


%% Plot M2 Peak over depth within the bathymetry

figure()

subplot(1,2,1)
contour(lon,lat,elv',-1000:-25:-3000,'-k','LineWidth',1,'Color',[0.5 0.5 0.5])
hold on
%plot(moorlon(5),moorlat(5),'o','LineWidth',2)
plot(moorlon(2),moorlat(2),'o','LineWidth',2)
plot(moorlon(1),moorlat(1),'o','LineWidth',2)
plot(moorlon(4),moorlat(4),'o','LineWidth',2)
plot(moorlon(3),moorlat(3),'o','LineWidth',2)
text(moorlon(1)+0.004,moorlat(1),moorname(1),'Linewidth',2)
text(moorlon(2)+0.004,moorlat(2),moorname(2),'Linewidth',2)
text(moorlon(3)-0.013,moorlat(3)+0.0015,moorname(3),'Linewidth',2)
text(moorlon(4)+0.004,moorlat(4)+0.001,moorname(4),'Linewidth',2)
%text(moorlon(5)+0.004,moorlat(5),moorname(5),'Linewidth',2)
title('Bathymetry and Mooring Locations')
xlims([-11.98 -11.8])
ylims([54.155 54.25])

subplot(1,2,2)
plot(mp2.M2peak,mp2.zAB)
hold on
plot(mavs1.M2peak,mavs1.zAB)
plot(tchain.M2peak,tchain.zAB)
plot(mp3.M2peak,mp3.zAB)
plot([0 0.5],[0 0],'k--')
legend('MP2','MAVS1','TCHAIN','MP3')
xlims([0 0.5])
title('Magnitude of M2 Peak')
xlabel('Power Spectral Density of M2 Peak [(m/s)^2/cpd]')
ylabel('Height Above Bottom [m]')



