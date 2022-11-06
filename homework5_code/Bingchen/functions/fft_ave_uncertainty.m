% The output of the function is averaged spectra from given segments number
% and their uncertainty 
% This function is used to load data(with constant data interval and without NaNs) 
%can only include 1 argument at a time. Note: data and time should have the
%same length


%There are three option for this function: 
% demean --> 1
% detrend -->2
% Hanning -->3

%Note: the default overlaping portion of data is 50%

function [spec_ave,uncertainty] = fft_ave_uncertainty(time,dt,data,sample_int,confi_int,option)


%% put data into segments



N = sample_int/2; 
nbin = floor((length(time)-N)/N);
for ii=1:nbin
    data_window{ii}= data(((ii-1)*N+1):((ii-1)*N+1)+2*N);
end 



%% 
for j = 1:length(data_window)

    if option == 1
        data_mod{j} = data_window{j}-mean(data_window{j});

        [f_mod(j,:),xs_mod(j,:),spec_mod(j,:)] = fft_data(data_mod{j},dt);

    elseif option == 2
        data_mod{j} = detrend(data_window{j});
        [f_mod(j,:),xs_mod(j,:),spec_mod(j,:)] = fft_data(data_mod{j},dt);


    elseif option == 3
        normf = sqrt(8/3);
        data_mod{j} = (hann(length(data_window{j}))*normf).*data_window{j};
        [f_mod(j,:),xs_mod(j,:),spec_mod(j,:)] = fft_data(data_mod{j},dt);
    else
        error('Please specify data modification option(must be either demean,detrend, or Hanning)')
    end 
end 


spec_ave = mean(spec_mod);


nu=2*nbin;
err_low = nu/chi2inv((confi_int/100)/2,nu);
err_high = nu/chi2inv((1-confi_int/100)/2,nu);
uncertainty(1,:) = err_high *spec_ave; % 1 raw is the lower bound 
uncertainty(2,:) = err_low *spec_ave; %2 raw is the upper bound 


end 








