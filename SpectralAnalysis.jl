module SpectralAnalysis
    using NetCDF, DataFrames
    using Statistics, LinearAlgebra, Distributions
    using DSP, Impute, FastTransforms

    function demean(data)
        μ = mean(data)
        return data .- μ
    end
    
    function trend(data, time)
        @assert length(time) == length(data)
        A = [ones(length(data)) time]
        b = A \ data
        return b
    end
    
    function detrend(data, time)
        b = trend(data, time)
        return data .- (b[1] .+ b[2] .* time)
    end
    
    """
        segment(df, var, ndays)
    
    df - the dataframe to segment
    var - the column of the dataframe to segment
    ndays - number of days per segment
    """
    function segment(data, time, T_segment; de_mean = false, de_trend = false, hanning = false, cosine = false)
        dt = time[2] - time[1] # the sample interval size [days]
        segment_size = iseven(round(T_segment / dt)) ? round(Int, T_segment / dt) : round(Int, T_segment / dt) - 1
    
        if de_trend
            @info "De-trending"
            data .= detrend(data, time)
        elseif de_mean
            @info "De-meaning"
            data .= demean(data)
        end
    
        segment_data = segment(data, segment_size)

        if hanning
            @info "Applying Hanning window"
            w = DSP.Windows.hanning(segment_size)
            segment_data = sqrt(8/3) .* w .* segment_data
        elseif cosine
            @info "Applying Cosine window"
            w = DSP.Windows.cosine(segment_size)
            segment_data = sqrt(2) .* w .* segment_data
        end
    
        return segment_data
    end
    
    # segment size is in # of points NOT days
    """
        segment(data, segment_size)
    
    Segments the data into 50% overlapping segments of length `segment_size`
    """
    function segment(data, segment_size::Int)
        @assert iseven(segment_size) # unclear how to overlap segments of odd size
        half_segment_size = Int(segment_size / 2)
        N = Int(floor(length(data) / half_segment_size)) # number of half-segments
        data = data[1:(N * half_segment_size)] # throw out some data at the end
        if iseven(N) # unshifted data captures last half segment
            data_unshift = data[1:(N * half_segment_size)]
            data_shift = data[(half_segment_size + 1):(end - half_segment_size)]
            @assert length(data_unshift) == length(data_shift) + segment_size
        else # shifted data captures last half segment
            data_unshift = data[1:(end - half_segment_size)]
            data_shift = data[(half_segment_size + 1):end]
            @assert length(data_unshift) == length(data_shift)
        end
    
        data_unshift = reshape(data_unshift, segment_size, Int(length(data_unshift) / segment_size))
        data_shift = reshape(data_shift, segment_size, Int(length(data_shift) / segment_size))
    
        data = zeros(segment_size, size(data_unshift)[2] + size(data_shift)[2])
        for i in 1:size(data_unshift)[2]
            data[:, 2*i - 1] .= data_unshift[:, i]
        end
    
        for i in 1:size(data_shift)[2]
            data[:, 2*i] .= data_shift[:, i]
        end
        return data
    end
    
    function frequencies(T, Δt)
        Δf = 1 / T # fundamental frequency f₁
        fₙ = 1 / (2Δt) # Nyquist frequency
        return Δf, fₙ
    end
    
    Φ(X, Δf) = norm.(X).^2 ./ Δf
    
    """
        spectral_uncertainty(n_segments; overlap = true)
    
    Returns the lower and upper 95% confidence interval coefficients for data with `n_segments` segments.
    If `overlap == true` then 50% overlap of segments is assumed.
    """
    function spectral_uncertainty(n_segments; overlap = true)
        ν = 2 * n_segments
        ν = overlap ? 0.9 * ν : ν # return 90% dof for 50% overlap
        chisq = Chisq(ν)
        lower = ν / quantile(chisq, 0.975)
        upper = ν / quantile(chisq, 0.025)
        
        return lower, upper
    end
    
    """
        spectrum(data, time; segment_size = 60, demean = true, detrend = true, hanning = false)
    
    data - data vector
    time - time vector
    sample_interval - number of days to include in a segment
    de_mean - boolean; if true, removes the mean from the data
    de_trend - boolean; if true, removes the trend from the data
    hanning - boolean; if true, applies a hanning window to segments
    """
    function spectrum(data, time, T_segment; de_mean = true, de_trend = true, hanning = false, cosine = false)
        # first, segment the data
        segment_data = segment(data, time, T_segment; de_mean = de_mean, de_trend = de_trend, hanning = hanning, cosine = cosine)
        segment_time = segment(time, time, T_segment; de_mean = false, de_trend = false, hanning =false, cosine = false)
        T = segment_time[end,1] - segment_time[1,1] # segment sample period
        Δt = time[2] - time[1] # the sample interval size [days]
        Δf, fₙ = frequencies(T, Δt)
        freq = Δf:Δf:(fₙ + Δf)
        spectra, uncertainty = spectrum(segment_data, Δf)
        return spectra, freq, uncertainty
    end
    
    function spectrum(segment_data, Δf)
        segment_size, n_segments = size(segment_data)
        X = fft(segment_data, 1) ./ segment_size

        # Check Parseval
        spectral_energy = mapslices(x -> mean(norm.(x).^2), X, dims = [1])
        data_energy = mapslices(x -> mean(x.^2) / segment_size, segment_data, dims = [1])
        @info "Parseval Check -- data energy: $(mean(data_energy)), spectral energy: $(mean(spectral_energy))"

        μ_X = mean(X, dims = 2)
        spectra = Φ(μ_X, Δf)[1:Int(segment_size / 2)]
        uncertainty = spectral_uncertainty(n_segments)
        return spectra, uncertainty
    end

end