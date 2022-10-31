% Lecture 10 coding exercise 


function [fft_amplitude,frequency_vector] = Okun_Lecture10(data,time,num_chunks)
    points_per_chunk = floor(length(data)/num_chunks); % Find nice number for points per chunk
    data_chunks = reshape(data,[points_per_chunk,num_chunks]);
    
    % For each chunk, calculate 
    for idx = 1:num_chunks
        fft_data(:,idx) = fft(data_chunks(:,idx));
        fft_amplitude(:,idx)=abs(fft_data(1:points_per_chunk/2+1,idx)).^2;
    end
    
    fft_amplitude = mean(fft_amplitude,2);
    time_span=(time(2)-time(1))*points_per_chunk; % Find time span in series
    df=1/time_span; % Frequency resolution
    fn=1/2/(time(2)-time(1));
    
    % Stolen from lecture 9
    frequency_vector=(0:df:fn)'; %frequency vector, cpd, goes from 0 to Nyquist.
    fft_amplitude = fft_amplitude / points_per_chunk.^2; % first correct for the MATLAB normalization
    fft_amplitude = fft_amplitude .* 2; %we threw out half of the spectrum; so correct for the lost variance.
    fft_amplitude = fft_amplitude / df; % this is then the definition of the spectrum

end