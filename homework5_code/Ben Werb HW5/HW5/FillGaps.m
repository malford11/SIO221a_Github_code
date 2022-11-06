function [df] = FillGaps(df)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Interpolate between points to fill gaps
%WSPD
idx = isnan(df.WSPD1);
df.WSPD1(idx) = interp1(df.t1(~idx),df.WSPD1(~idx),df.t1(idx),'linear');
idx = isnan(df.WSPD2);
df.WSPD2(idx) = interp1(df.t2(~idx),df.WSPD2(~idx),df.t2(idx),'linear');
idx = isnan(df.WSPD3);
df.WSPD3(idx) = interp1(df.t3(~idx),df.WSPD3(~idx),df.t3(idx),'linear');
%VWND
idx = isnan(df.VWND1);
df.VWND1(idx) = interp1(df.t1(~idx),df.VWND1(~idx),df.t1(idx),'linear');
idx = isnan(df.VWND2);
df.VWND2(idx) = interp1(df.t2(~idx),df.VWND2(~idx),df.t2(idx),'linear');
idx = isnan(df.WSPD3);
df.VWND3(idx) = interp1(df.t3(~idx),df.VWND3(~idx),df.t3(idx),'linear');
%UWND
idx = isnan(df.UWND1);
df.UWND1(idx) = interp1(df.t1(~idx),df.UWND1(~idx),df.t1(idx),'linear');
idx = isnan(df.UWND2);
df.UWND2(idx) = interp1(df.t2(~idx),df.UWND2(~idx),df.t2(idx),'linear');
idx = isnan(df.UWND3);
df.UWND3(idx) = interp1(df.t3(~idx),df.UWND3(~idx),df.t3(idx),'linear');
%remove nans at end of time 1
idx = find(isnan(df.WSPD1));
df.WSPD1(idx) = []; %drop extra NaNs at end of record
idx = find(isnan(df.VWND1));
df.VWND1(idx) = [];
idx = find(isnan(df.UWND1));
df.UWND1(idx) = [];
df.t1(idx) = [];
%Remake df.WSPD
df.WSPD = [df.WSPD1 , df.WSPD2 , df.WSPD3];
%Remake df.VWND
df.VWND = [df.VWND1 , df.VWND2 , df.VWND3];
%Remake df.UWND
df.UWND = [df.UWND1 , df.UWND2 , df.UWND3];
%Remake time
df.t = [df.t1 ; df.t2 ; df.t3];
df.t_dnum=df.t+datenum(1950,1,1); %days since 1950
end

