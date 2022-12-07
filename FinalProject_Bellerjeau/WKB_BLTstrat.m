function [mooring] = WKB_BLTstrat(mooring)
%%%%% This function adds stretched and scaled velocity and z coordinates to
%%%%% your pre-existing structure 'mooring'. Run mooring_load first to
%%%%% convert your.nc data file into a structure this function will work
%%%%% with. Stratification is currently set to the BLT Canyon.
%%%%% 
%%%%% Inputs: struct 'mooring'
%%%%% Returns: struct 'mooring' with 4 new fields added 
%%%%%       - mooring.uSc = scaled U velocity
%%%%%       - mooring.vSc = scaled V velocity
%%%%%       - mooring.wSc = scaled W velocity
%%%%%       - mooring.zSc = scaled z coordinates
%%%%%       - mooring.N2 = N^2 profile at mooring depths
%%%%%
%%%%% by Charlotte Bellerjeau 10/3/22

%load CTD transect data (T4) for N^2 profile
load('n2_ctds.mat')
load('z_ctds.mat')
load('CTDcrosscanyon.mat')
load('MP1.mat')

[m,n] = size(mooring.u);

%% Obtain N^2 profile for best fit

%averaging CTD data from T4 transect
Zout = zeros(length(n2_ctds(:,1)),1);
N2out = zeros(length(n2_ctds(:,1)),1);
for i = 1:length(z_ctds(:,1))
    N2 = n2_ctds(i,:);
    N2(isnan(N2)) = [];
    N2out(i) = mean(N2);
    Z = z_ctds(i,:);
    Z(isnan(Z)) = [];
    Zout(i) = mean(Z);
end

%trim mixed layer for nice profile
Zout = Zout(300:3000);
N2out = N2out(300:3000);

%smoothing
N2smooth = smoothdata(N2out,'movmean',100);

%plot N^2 profile to check smoothing
% figure()
% plot(N2out,Zout)
% hold on
% plot(N2smooth,Zout)
% plot(MPall.N2mean,MPall.z)
% plot(CTD.n2m,CTD.z)
% hold off
% axis ij
% title('N^2 Profile')
% legend('Averaged T4 CTDs','Smoothed N2 Profile','Mooring', 'Cross Canyon FCTD')
% xlabel('N^2 [1/s]')
% ylabel('Depth [m]')

%% Scaling Velocity

%interpolate to depths at which we have data points
N2 = (interp1(Zout,N2smooth,mooring.z));
N2(find(isnan(N2)))=N2(find(isnan(N2))+1);     %remove NaNs

N = sqrt(N2);
Nbar = mean(N);

%scale u and v by stratification
stratscale = (N./Nbar).^(1/2);
for i = 1:n
    mooring.uSc(:,i) = mooring.u(:,i)./ stratscale(i);
    mooring.vSc(:,i) = mooring.v(:,i)./stratscale(i);
    mooring.wSc(:,i) = mooring.w(:,i)./stratscale(i);
    mooring.upSc(:,i) = mooring.up(:,i)./stratscale(i);
    mooring.crSc(:,i) = mooring.cr(:,i)./stratscale(i);
end

%% Stretching Vertical Coordinate
mooring.zSc = cumtrapz(mooring.z,N./Nbar) + mooring.z(1);

%plot stretched and scaled velocity alongside original velocity
% figure()
% plot(mooring.u(1,:),mooring.z)
% hold on
% plot(mooring.uSc(1,:),mooring.zSc)
% hold off
% axis ij
% title('WKB Stretch and Scaling for Stratification')
% xlabel('U Velocity')
% ylabel('Depth [m]')
% legend('Unscaled','Scaled')

%% Add N^2 to mooring structure
mooring.N2 = N2;

end
