using NetCDF, DataFrames
using Plots
using Statistics, LinearAlgebra, Distributions
using DSP, Impute, FastTransforms

FILE_PATH = joinpath(@__DIR__, "HW5_data/")
FIGURE_PATH = joinpath(@__DIR__, "hw5_figures/")

########################### Function Definitions #####################################
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
function segment(data, time, ndays; de_mean = false, de_trend = false, hanning = false)
    dt = time[2] - time[1] # the sample interval size [days]
    segment_size = iseven(round(ndays / dt)) ? round(Int, ndays / dt) : round(Int, ndays / dt) - 1
    segment_data = segment(data, segment_size)
    segment_time = segment(time, segment_size)

    if de_trend
        @info "De-trending"
        segment_data = map(zip(eachcol(segment_data), eachcol(segment_time))) do (x, t)
            detrend(x, t)
        end
        segment_data = hcat(Tuple(segment_data)...)
    elseif de_mean
        @info "De-meaning"
        segment_data = mapslices(demean, segment_data, dims = 1)
    end

    if hanning
        @info "Applying Hanning window"
        w = DSP.Windows.hanning(segment_size)
        segment_data = sqrt(8/3) .* w .* segment_data
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
    ν = overlap ? 0.9 * ν : ν # return 90% dof for hanning overlap
    chisq = Chisq(ν)
    lower = ν / quantile(chisq, 0.025)
    upper = ν / quantile(chisq, 0.975)
    
    return lower, upper
end

"""
    spectrum(data, time; sample_interval = 60, demean = true, detrend = true, hanning = false)

data - data vector
time - time vector
sample_interval - number of days to include in a segment
de_mean - boolean; if true, removes the mean from the data
de_trend - boolean; if true, removes the trend from the data
hanning - boolean; if true, applies a hanning window to segments
"""
function spectrum(data, time; segment_size_days = 60, de_mean = true, de_trend = true, hanning = false)
    # first, segment the data
    segment_data = segment(data, time, segment_size_days; de_mean = de_mean, de_trend = de_trend, hanning = hanning)
    segment_time = segment(time, time, segment_size_days; de_mean = false, de_trend = false, hanning =false)
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
    μ_X = mean(X, dims = 2)
    spectral_energy = mapslices(x -> sum(norm.(x).^2), X, dims = [1])
    data_energy = mapslices(x -> mean(x.^2), segment_data, dims = [1])
    @info "Parseval Check -- data energy: $(mean(data_energy)), spectral energy: $(mean(spectral_energy))"
    spectra = Φ(μ_X, Δf)[1:Int(segment_size / 2)]
    uncertainty = spectral_uncertainty(n_segments)
    return spectra, uncertainty
end

#########################################################################

# 1. Load data & Preliminary analysis

file1 = joinpath(FILE_PATH, "OS_T8S110W_DM134A-20150425_D_WIND_10min.nc")
file2 = joinpath(FILE_PATH, "OS_T8S110W_DM183A-20160321_D_WIND_10min.nc")
file3 = joinpath(FILE_PATH, "OS_T8S110W_DM231A-20170606_D_WIND_10min.nc")

# read in time data [days since 1950-01-01T00:00:00Z]
time1 = ncread(file1, "TIME")
time2 = ncread(file2, "TIME")
time3 = ncread(file3, "TIME")
time = vcat(time1, time2, time3)

@info "Time gaps: $(24 * (time2[1] - time1[end])) hrs, $(24 * (time3[1] - time2[end])) hrs"

plot(time[2:end], diff(time),
     title = "Sample time intervals",
     xlabel = "Days since 1950-01-01T00:00:00Z",
     ylabel = "Sampling Gap [Days]",
     legend = false,
    )
# savefig(joinpath(FIGURE_PATH, "diff.png"))

# wind speed [m/s]
wspd1 = ncread(file1, "WSPD")
wspd2 = ncread(file2, "WSPD")
wspd3 = ncread(file3, "WSPD")
wspd = vec(hcat(wspd1, wspd2, wspd3))

plot(time, wspd,
     title = "Raw Wind Speed Data",
     xlabel = "Days since 1950-01-01T00:00:00Z",
     ylabel = "Wind speed [m/s]",
     legend = false,
    )
# savefig(joinpath(FIGURE_PATH, "wspd_raw.png"))

# zonal wind [m/s]
uwnd1 = ncread(file1, "UWND")
uwnd2 = ncread(file2, "UWND")
uwnd3 = ncread(file3, "UWND")
uwnd = vec(hcat(uwnd1, uwnd2, uwnd3))

plot(time, uwnd,
     title = "Raw Zonal Wind Data",
     xlabel = "Days since 1950-01-01T00:00:00Z",
     ylabel = "Wind speed [m/s]",
     legend = false,
    )
# savefig(joinpath(FIGURE_PATH, "uwnd_raw.png"))

# meridional wind [m/s]
vwnd1 = ncread(file1, "VWND")
vwnd2 = ncread(file2, "VWND")
vwnd3 = ncread(file3, "VWND")
vwnd = vec(hcat(vwnd1, vwnd2, vwnd3))

plot(time, vwnd,
     title = "Raw Meridional Wind Data",
     xlabel = "Days since 1950-01-01T00:00:00Z",
     ylabel = "Wind speed [m/s]",
     legend = false,
    )
# savefig(joinpath(FIGURE_PATH, "vwnd_raw.png"))

df = DataFrame(
    :time => time,
    :wspd => wspd,
    :uwnd => uwnd,
    :vwnd => vwnd
)

# impute missing data
df = Impute.declaremissings(df, values = (-999.0f0))
Impute.interp!(df)
disallowmissing!(df)

df_2015 = view(df, 1:length(time1), :)
df_2016 = view(df, (1:length(time2)) .+ length(time1), :)
df_2017 = view(df, (1:length(time3)) .+ (length(time1) + length(time2)), :)

plot(df.time, df.wspd,
     title = "Filled Windspeed Data",
     xlabel = "Days since 1950-01-01T00:00:00Z",
     ylabel = "Wind speed [m/s]",
     legend = false,
    )
# savefig(joinpath(FIGURE_PATH, "wspd_fill.png"))

plot(df.time, df.uwnd,
     title = "Filled Zonal Wind Data",
     xlabel = "Days since 1950-01-01T00:00:00Z",
     ylabel = "Wind speed [m/s]",
     legend = false,
   )
# savefig(joinpath(FIGURE_PATH, "uwnd_fill.png"))

plot(df.time, df.vwnd,
     title = "Filled Meridional Wind Data",
     xlabel = "Days since 1950-01-01T00:00:00Z",
     ylabel = "Wind speed [m/s]",
     legend = false,
    )
# savefig(joinpath(FIGURE_PATH, "vwnd_fill.png"))

# 2. segment data
cols = NamedTuple(zip(Symbol.(names(df)), Symbol.(names(df))))
segment_2015 = map(cols) do col
    segment(getproperty(df_2015, col), df_2015.time, 60)
end
pt_per_seg_2015, N_seg_2015 = size(segment_2015.time)
Δf_2015, fₙ_2015 = frequencies(segment_2015.time[end, 1] - segment_2015.time[1, 1], segment_2015.time[2] - segment_2015.time[1])
@info "____________2015____________"
@info "n_segments = $(N_seg_2015), pt_per_segment = $(pt_per_seg_2015)"
@info "f₁ = $(Δf_2015), fₙ = $(fₙ_2015)"

segment_2016 = map(cols) do col
    segment(getproperty(df_2016, col), df_2016.time, 60)
end
pt_per_seg_2016, N_seg_2016 = size(segment_2016.time)
Δf_2016, fₙ_2016 = frequencies(segment_2016.time[end, 1] - segment_2016.time[1, 1], segment_2016.time[2] - segment_2016.time[1])
@info "____________2016____________"
@info "n_segments = $(N_seg_2016), pt_per_segment = $(pt_per_seg_2016)"
@info "f₁ = $(Δf_2016), fₙ = $(fₙ_2016)"

segment_2017 = map(cols) do col
    segment(getproperty(df_2017, col), df_2017.time, 60)
end
pt_per_seg_2017, N_seg_2017 = size(segment_2017.time)
Δf_2017, fₙ_2017 = frequencies(segment_2017.time[end, 1] - segment_2017.time[1, 1], segment_2017.time[2] - segment_2017.time[1])
@info "____________2017____________"
@info "n_segments = $(N_seg_2017), pt_per_segment = $(pt_per_seg_2017)"
@info "f₁ = $(Δf_2017), fₙ = $(fₙ_2017)"

# 3. raw vs detrend vs detrend + Hanning

# raw
Φ_demean, frange, (lower_2015, upper_2015) = spectrum(df_2015.wspd, df_2015.time; de_mean = true, de_trend = false)
a, b = 1, length(frange)
plot(frange[a:b], Φ_demean[a:b],
     title = "Wind Speed Spectrum",
     ylabel = "Φₓ (m/s)²/cpd)",
     xlabel = "frequency [cycles/day]",
     labels = "De-meaned",
     xaxis=:log,
     yaxis=:log,
     legend = :bottomright
    )

# detrend
Φ_detrend, _, _ = spectrum(df_2015.wspd, df_2015.time)
plot!(frange[a:b], Φ_detrend[a:b], labels = "De-trended",)

# detrend + hanning window
Φ_hanning, _, _ = spectrum(df_2015.wspd, df_2015.time; hanning = true)
plot!(frange[a:b], Φ_hanning[a:b], labels = "De-trended + Hanning")

# spectral uncertainties
plot(frange[a:b], Φ_demean[a:b],
     ribbon = (Φ_demean[a:b] .- (Φ_demean[a:b] .* lower_2015), (Φ_demean[a:b] .* upper_2015) .- Φ_demean[a:b]),
     labels = "De-meaned",
     seriescolor = "green",
     title = "2015 Wind Speed Spectrum",
     ylabel = "Φₓ (m/s)²/cpd)",
     xlabel = "frequency [cycles/day]",
     xaxis=:log,
     yaxis=:log,
     legend = :bottomleft,
     yrange = (10^-10, 10)
    )
plot!(frange[a:b], Φ_detrend[a:b],
     ribbon = (Φ_detrend[a:b] .- (Φ_detrend[a:b] .* lower_2015), (Φ_detrend[a:b] .* upper_2015) .- Φ_detrend[a:b]),
     labels = "De-trended",
     seriescolor = "red",
    )
plot!(frange[a:b], Φ_hanning[a:b],
    ribbon = (Φ_hanning[a:b] .- (Φ_hanning[a:b] .* lower_2015), (Φ_hanning[a:b] .* upper_2015) .- Φ_hanning[a:b]),
    labels = "De-trended + Hanning",
    seriescolor = "blue",
   )

# 4. multi-year record vs single year record
segment_time = segment(df_2015.time, df_2015.time, 60; de_mean = false, de_trend = false, hanning = false)
wspd_segment_2015 = segment(df_2015.wspd, df_2015.time, 60; de_trend = true, hanning = true)
wspd_segment_2016 = segment(df_2016.wspd, df_2016.time, 60; de_trend = true, hanning = true)
wspd_segment_2017 = segment(df_2017.wspd, df_2017.time, 60; de_trend = true, hanning = true)
full_wspd = hcat(wspd_segment_2015, wspd_segment_2016, wspd_segment_2017)

uwnd_segment_2015 = segment(df_2015.uwnd, df_2015.time, 60; de_trend = true, hanning = true)
uwnd_segment_2016 = segment(df_2016.uwnd, df_2016.time, 60; de_trend = true, hanning = true)
uwnd_segment_2017 = segment(df_2017.uwnd, df_2017.time, 60; de_trend = true, hanning = true)
full_uwnd = hcat(uwnd_segment_2015, uwnd_segment_2016, uwnd_segment_2017)

vwnd_segment_2015 = segment(df_2015.vwnd, df_2015.time, 60; de_trend = true, hanning = true)
vwnd_segment_2016 = segment(df_2016.vwnd, df_2016.time, 60; de_trend = true, hanning = true)
vwnd_segment_2017 = segment(df_2017.vwnd, df_2017.time, 60; de_trend = true, hanning = true)
full_vwnd = hcat(vwnd_segment_2015, vwnd_segment_2016, wspd_segment_2017)

T = segment_time[end,1] - segment_time[1,1] # segment sample period
Δt = segment_time[2] - segment_time[1] # the sample interval size [days]
Δf, fₙ = frequencies(T, Δt)
freq = Δf:Δf:(fₙ + Δf)

Φ_wspd, (lower_all_yr, upper_all_yr) = spectrum(full_wspd, Δf)
Φ_uwnd, _ = spectrum(full_uwnd, Δf)
Φ_vwnd, _ = spectrum(full_vwnd, Δf)

# 5. 2015 vs all years - windspeed
a1, b1 = 1, length(frange)
plot(frange[a1:b1], Φ_hanning[a1:b1],
    title = "Wind Speed Spectra",
    ylabel = "Φₓ (m/s)²/cpd)",
    xlabel = "frequency [cycles/day]",
    ribbon = (Φ_hanning[a1:b1] .- (Φ_hanning[a1:b1] .* lower_2015), (Φ_hanning[a1:b1] .* upper_2015) .- Φ_hanning[a1:b1]),
    labels = "2015",
    seriescolor = "blue",
    xaxis=:log,
    yaxis=:log,
    yrange = (10^-10, 10),
    xrange = (0.5, 5)
   )
a2, b2 = 1, length(frange_all_yr)
plot!(frange_all_yr[a2:b2], Φ_wspd[a2:b2],
    ribbon = (Φ_wspd[a2:b2] .- (Φ_wspd[a2:b2] .* lower_all_yr), (Φ_wspd[a2:b2] .* upper_all_yr) .- Φ_wspd[a2:b2]),
    labels = "2015-2017",
    seriescolor = "red",
  )

# 6. all wind data spectra for all years
plot(frange_all_yr[a2:b2], Φ_wspd[a2:b2],
    ribbon = (Φ_wspd[a2:b2] .- (Φ_wspd[a2:b2] .* lower_all_yr), (Φ_wspd[a2:b2] .* upper_all_yr) .- Φ_wspd[a2:b2]),
    title = "Wind Spectra\n2015-2017",
    ylabel = "Φₓ (m/s)²/cpd)",
    xlabel = "frequency [cycles/day]",
    labels = "wind speed",
    legend = :bottomleft,
    seriescolor = "green",
    xaxis=:log,
    yaxis=:log,
    yrange = (10^-7, 10),
    xrange = (0.5, 5)
  )

plot!(frange_all_yr[a2:b2], Φ_uwnd[a2:b2],
  ribbon = (Φ_uwnd[a2:b2] .- (Φ_uwnd[a2:b2] .* lower_all_yr), (Φ_uwnd[a2:b2] .* upper_all_yr) .- Φ_uwnd[a2:b2]),
#   title = "Zonal Wind Spectrum\n2015-2017",
  ylabel = "Φₓ (m/s)²/cpd)",
  xlabel = "frequency [cycles/day]",
  labels = "zonal wind",
  seriescolor = "blue",
  xaxis=:log,
  yaxis=:log,
)

plot!(frange_all_yr[a2:b2], Φ_vwnd[a2:b2],
    ribbon = (Φ_vwnd[a2:b2] .- (Φ_vwnd[a2:b2] .* lower_all_yr), (Φ_vwnd[a2:b2] .* upper_all_yr) .- Φ_vwnd[a2:b2]),
    # title = "Meridional Wind Spectrum\n2015-2017",
    ylabel = "Φₓ (m/s)²/cpd)",
    xlabel = "frequency [cycles/day]",
    labels = "meridional wind",
    seriescolor = "red",
    xaxis=:log,
    yaxis=:log,
  )