function [df] = Structify(data1,data2,data3)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
df = struct;
df.t1 = ncread(data1,'TIME');
df.t2 = ncread(data2,'TIME');
df.t3 = ncread(data3,'TIME');
df.t = [df.t1 ; df.t2 ; df.t3];

df.t_dnum=df.t+datenum(1950,1,1); %days since 1950
df.t1_dnum=df.t1+datenum(1950,1,1);
df.t2_dnum=df.t2+datenum(1950,1,1);
df.t3_dnum=df.t3+datenum(1950,1,1);

df.H1 = ncread(data1,'HEIGHT');
df.H2 = ncread(data2,'HEIGHT');
df.H3 = ncread(data3,'HEIGHT');
df.H = [df.H1 ; df.H2 ; df.H3];

df.WSPD1 = ncread(data1,'WSPD');
df.WSPD2 = ncread(data2,'WSPD');
df.WSPD3 = ncread(data3,'WSPD');
df.WSPD = [df.WSPD1 , df.WSPD2 , df.WSPD3];

df.VWND1 = ncread(data1,'VWND');
df.VWND2 = ncread(data2,'VWND');
df.VWND3 = ncread(data3,'VWND');
df.VWND = [df.VWND1 , df.VWND2 , df.VWND3];

df.UWND1 = ncread(data1,'UWND');
df.UWND2 = ncread(data2,'UWND');
df.UWND3 = ncread(data3,'UWND');
df.UWND = [df.VWND1 , df.VWND2 , df.VWND3];

df.names = fieldnames(df);
df.names = df.names([16 20 24]);
df.units = {'m/s', 'm/s', 'm/s'};
end

