% function for homework 7, written by Bingchen Liu
% This code first put data in to chunks and calculate the fft for each
% chunk and take the average 

% note that the data is detrended 

function [spec_ave,spec,f,data] = Liu_lecture10_func(time,data_input,sample_int)

dt=mean(diff(time),'omitnan');

n_int = floor(length(data_input)/sample_int);
data_mod = data_input(1:n_int*sample_int);

data = reshape(data_mod,[n_int,sample_int]);

% f = zeros(n_int,length(data(1,:)));
% spec =  zeros(n_int,length(data(1,:)));
    for j = 1:n_int

       [f(j,:),spec(j,:)] = fft_data(detrend(data(j,:)),dt);
    end 
    
spec_ave = mean(spec,1);

end 