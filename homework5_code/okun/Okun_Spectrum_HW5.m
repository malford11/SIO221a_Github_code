% Okun function for HW5

function [fft_amplitude,frequency_vector] = Okun_Spectrum_HW5(data,time,detrend_demean,window)

    if ~exist('window','var')
        window = 0;
    end

    if ~exist('detrend','var')
        detrend_demean = 0;
    end

    data = fillmissing(data,'linear');
    if detrend_demean
        data=detrend(data);
    end

    if window
        days_per_segment = 60;
        points_per_chunk = days_per_segment*26*6;
        hann_window = 8/3*hann(points_per_chunk);
        num_chunks = floor(length(data)/(points_per_chunk/2));
        for idx = 1:(num_chunks-1)
            data_chunks(:,idx) = data((1+(idx-1)*points_per_chunk/2):(points_per_chunk+(idx-1)*points_per_chunk/2));
            % Apply hanning window
            data_chunks(:,idx) = detrend(data_chunks(:,idx));
            data_chunks(:,idx) = data_chunks(:,idx).*hann_window;
        end



        % For each chunk, calculate
        for idx = 1:(num_chunks-1)
            fft_data(:,idx) = fft(data_chunks(:,idx));
            fft_amplitude(:,idx)=abs(fft_data(1:points_per_chunk/2+1,idx)).^2;
        end
    else
        points_per_chunk = length(data);
        data_chunks = data;
        fft_data = fft(data_chunks);
        fft_amplitude=abs(fft_data(1:points_per_chunk/2+1)).^2;
    end

    if window
        fft_amplitude = mean(fft_amplitude,2);
    end

    dt=mean(diff(time),'omitnan');
    fft_amplitude = mean(fft_amplitude,2);
    time_span=dt*points_per_chunk; % Find time span in series
    df=1/time_span; % Frequency resolution
    fn=1/2/dt;
    

    frequency_vector=(0:df:fn)'; %frequency vector, cpd, goes from 0 to Nyquist.
    fft_amplitude = fft_amplitude / points_per_chunk.^2; % first correct for the MATLAB normalization
    fft_amplitude = fft_amplitude .* 2; %we threw out half of the spectrum; so correct for the lost variance.
    fft_amplitude = fft_amplitude / df; % this is then the definition of the spectrum
    
    figure
    loglog(frequency_vector,fft_amplitude)
    if window
        title('Power Spectral Density Windowed')
    else
        title('Power Spectral Density')
    end
    xlabel('Frequency (cpd)')
    ylabel('\Phi_{Data}  ((units)^2 / cpd)')
end