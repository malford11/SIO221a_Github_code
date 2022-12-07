function [mooring] = mooring_load(file, trim)

%%%%%% Run this script to load BLT Mooring MP3 dataset and trim NaNNed rows
%%%%%% and columns. Also creates useful colormaps, zeros out remaining
%%%%%% NaNs, and assigns the size of the trimmed dataset to variables


%% Load Data and Trim NaNs

mooring=ncinfo(file);
%mooring.z=ncread(file,'depth');  %some mooring files have 'depth' while others have 'z'
mooring.u=ncread(file,'u');
mooring.v=ncread(file,'v');
mooring.w=ncread(file,'w');
mooring.t=ncread(file,'temperature');       %temperature
mooring.p=ncread(file,'pressure');          %pressure
mooring.time=ncread(file,'time');
mooring.time=cast(mooring.time,'double');
% % MP3: 
%mooring.time = mooring.time./1e6./60./60./24;  % time is now in days since 2021-10-15
% MP2 and 1:
%mooring.time = mooring.time./60./24;  % time is now in days since 2021-07-07 10:15am and 2021-06-28 15:15:00'

% Trimming NaNs
% trim = [1 10 27248; 1 2 35]

Mats = {'u','v','w'};
for i = 1:length(Mats)
    mooring.(Mats{i})(1:trim(1,2),:) =[];
    mooring.(Mats{i})(trim(1,3):end,:)=[];
    mooring.(Mats{i})(:,1:trim(2,2))=[];
    mooring.(Mats{i})(:,trim(2,3):end)=[];
end
Scals = {'time','t','p'};
for i = 1:length(Scals)
    mooring.(Scals{i})(1:trim(1,2))=[];
    mooring.(Scals{i})(trim(1,3):end)=[];
end
% mooring.z(1:trim(2,2))=[];
% mooring.z(trim(2,3):end)=[];

%zero out remaining NaNs for now, may interpolate later
for i=1:length(Mats)
    mooring.(Mats{i})((isnan(mooring.(Mats{i}))))=0;
end

%rotate coordinates to be up-canyon/cross-canyon
[mooring.up,mooring.cr]=RotateCCW(mooring.u,mooring.v,7*pi/4);

