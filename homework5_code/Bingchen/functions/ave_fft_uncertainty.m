%Note: this function is used in HW5 specifically, not the final version of
%problem 7. For the function of problem 7, the name is fft_ave_uncertainty
% Input of this function is segmented data (with same length)
%confi_int(confidence interval) is in percentage

function [f,spec_ave,uncertainty] = ave_fft_uncertainty(data_seg,dt,confi_int)

%dt = 0.006944444441615; % in days
for j = 1:length(data_seg)

    %detrended data 
    data_detre{j} = detrend(data_seg{j});
    [f_temp.detre(j,:),xs.detre(j,:),spec.detre(j,:)] = fft_data(data_detre{j},dt);

    %Hanning window applied 
    normf = sqrt(8/3);
    data_hanning{j} = (hann(length(data_seg{j}))*normf).*detrend(data_seg{j});
    [f_temp.hanning(j,:),xs.hanning(j,:),spec.hanning(j,:)] = fft_data(data_hanning{j},dt);
end 

f.detre = f_temp.detre(1,:);
f.hanning = f_temp.hanning(1,:);

spec_ave.detre = mean(spec.detre);
spec_ave.hanning = mean(spec.hanning);


nu=2*length(data_seg);
err_low = nu/chi2inv((confi_int/100)/2,nu);
err_high = nu/chi2inv((1-confi_int/100)/2,nu);
uncertainty.hanning(1,:) = err_high *spec_ave.hanning; % 1 raw is the lower bound 
uncertainty.hanning(2,:) = err_low *spec_ave.hanning; %2 raw is the upper bound 

uncertainty.detre(1,:) = err_high *spec_ave.detre; % 1 raw is the lower bound 
uncertainty.detre(2,:) = err_low *spec_ave.detre; %2 raw is the upper bound 


end 
