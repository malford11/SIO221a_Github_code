function [a] = InterpolateData(data, time)
%function for removing NaN values from data and using interp1 to replace
%them
%takes data and time vector as input
%outputs interpolated data

ind = isnan(data);
data(ind) = [];
t = time;
t(ind) = [];
a = interp1(t, data, time);