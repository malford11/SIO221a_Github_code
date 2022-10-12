function [df] = Scripps_Pier(df)
%Scripps Pier grabs time, temp, and pressure and puts them in a struct.
%Data files from https://sccoos.org/autoss/

time = ncread("scripps-pier-automated-shore-sta_1dd8_b354_6796.nc","time");
t = ncread("scripps-pier-automated-shore-sta_1dd8_b354_6796.nc", "sea_water_temperature");
p = ncread("scripps-pier-automated-shore-sta_f6dc_49c5_7e47.nc",'sea_water_pressure');
df = struct;
df.time = time;
df.t = t;
date0=datenum(1970,1,1); %put in correct time format
dnum=double(df.time)/3600/24+date0;
df.time = dnum;
df.p = p;
idx1 = datenum(2005,1,1);
idx2 = datenum(2023,1,1);
idx = find(df.time >= idx1 & df.time <= idx2);
df.time = df.time(idx);
df.t = df.t(idx);
df.p = df.p(idx);
df.datetime = datevec(df.time); %easy to view and sort dates

%remove outliers:
outliersup = 3.5 * abs(std(df.t));
outliersdown = -3.5 * abs(std(df.t));
idx = find(df.t >= outliersup & df.t <= outliersdown);
df.t(idx) = [];

outliersup = 3.5 * abs(std(df.p));
outliersdown = -3.5 * abs(std(df.p));
idx = find(df.p >= outliersup & df.p <= outliersdown);
df.p(idx) = [];

end

