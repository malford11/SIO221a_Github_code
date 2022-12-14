% Okun Lecture 10 coding exercise
% 11/9/2022
% Creates a spectrum and frequency vector given data, time, and number of chunks.
% Output is is fft coefficients with units (data_units)^2/(1/time_units)
% and frequency vector with units 1/time_units

function [fft_amplitude,frequency_vector] = Okun_Lecture10_Hanning(data,time,num_chunks)
    
    % Fill gaps in data and detrend    
    data = fillmissing(data,'linear');
    data=detrend(data);
    
    % Find nice number for points per chunk
    points_per_chunk = floor(length(data)/num_chunks);
    
    % Make length of chunks even if not
    if mod(points_per_chunk,2)
        points_per_chunk = points_per_chunk - 1;
    end
    
    % Reshape data into chunks
    data_chunks = reshape(data(1:points_per_chunk*num_chunks),[],num_chunks);
    
    % Remove mean from data
    data_chunks = data_chunks - mean(data_chunks);

    % For each chunk, calculate 
    for idx = 1:num_chunks
        fft_data(:,idx) = fft(data_chunks(:,idx));
        fft_amplitude(:,idx)=abs(fft_data(1:points_per_chunk/2+1,idx)).^2;
    end
    
    % Create frequency vector
    dt=mean(diff(time),'omitnan');
    fft_amplitude = mean(fft_amplitude,2);
    time_span=dt*points_per_chunk; % Find time span in series
    df=1/time_span; % Frequency resolution
    fn=1/2/dt;
    
    % Stolen from lecture 9
    frequency_vector=(0:df:fn)'; %frequency vector, cpd, goes from 0 to Nyquist.
    fft_amplitude = fft_amplitude / points_per_chunk.^2; % first correct for the MATLAB normalization
    fft_amplitude = fft_amplitude .* 2; %we threw out half of the spectrum; so correct for the lost variance.
    fft_amplitude = fft_amplitude / df; % this is then the definition of the spectrum
    
    variance = mean(data.^2,"omitnan");
    sum_spec = sum(fft_amplitude)*df;
    ratio = sum_spec / variance;

    disp(['The variance of the data is ' num2str(variance) ' and the variance of the spectrum is ' num2str(sum_spec)])
    disp(['The ratio of the variance of the data to the spectrum is ' num2str(ratio)])
end