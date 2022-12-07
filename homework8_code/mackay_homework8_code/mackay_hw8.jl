using Impute, Statistics, LinearAlgebra
using Plots
using NetCDF, DataFrames

include("SpectralAnalysis.jl")

FILE_PATH  = joinpath(@__DIR__, "HW5_data/")
FIGURE_PATH = joinpath(@__DIR__, "HW8_figures")

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

# Normalization & Parseval Check
dt = df_2015.time[2] - df_2015.time[1]
N = length(df_2015.time)
Φ, freq, (lower, upper) = SpectralAnalysis.spectrum(df_2015.wspd, df_2015.time, N * dt)

# Check Parseval
@info var(df_2015.wspd) / (sum(Φ) * freq[1])

Φ_hann, freq, (lower, upper) = SpectralAnalysis.spectrum(df_2015.wspd, df_2015.time, N / 15 * dt, hanning = true)
@info var(df_2015.wspd) / (sum(Φ_hann) * freq[1])

function f_alias(f_nyquist, f_signal)
    M = floor(f_signal / f_nyquist)
    isodd(M) ? (return ((M+1) * f_nyquist - f_signal)) : (return (f_signal - M * f_nyquist))
end
 
# Spectra of aliased signals
periods = (S1 = 24.00 / 24., TwoN2 = 12.9054 / 24., N2 = 12.6583 / 24., M2 = 12.4206 / 24., S2 = 12.00 / 24., K2 = 11.9672/ 24.) # [hrs]
sample_periods = (fast = .99349, science = 20.86455) # [cpd]
for (f, f_name) in zip(sample_periods, keys(sample_periods))
    @info "________$(f_name)_______"
    freqs = zeros(length(periods))
    alias_freqs = similar(freqs)
    fn = 1 / (2 * f)
    @info "Nyquist Freq: $(fn) cpd"
    for i in 1:length(periods)
        T, T_name = periods[i], keys(periods)[i]
        freqs[i] = 1 / T
        alias_freq = f_alias(fn, freqs[i])
        alias_freqs[i] = alias_freq
        @info "$T_name: $(alias_freq)"
    end
    sort!(freqs)
    sort!(alias_freqs)
    @info "Minimum frequency separation w/o alias: $(minimum(diff(freqs))) cpd"
    @info "Necessary operation duration w/o alias: $(1 / minimum(diff(freqs))) days"
    @info "Minimum frequency separation w/ alias: $(minimum(diff(alias_freqs))) cpd"
    @info "Necessary operation duration w/ alias: $(1 / minimum(diff(alias_freqs))) days"
end

# Spectra of aliased signals
dt = df_2015.time[2] - df_2015.time[1]
T_segment = 60
Φ, freq, (lower, upper) = SpectralAnalysis.spectrum(df_2015.wspd, df_2015.time, T_segment; hanning = true)

time_subsample = df_2015.time[1:40:length(df_2015.time)]
wspd_subsample = df_2015.wspd[1:40:length(df_2015.time)]
Φ_subsample, freq_subsample, (lower_subsample, upper_subsample) = SpectralAnalysis.spectrum(wspd_subsample, time_subsample, T_segment; hanning = true)

plot(freq, Φ,
    xaxis = :log,
    yaxis = :log,
    legend = :bottomright,
    labels = "standard spectrum",
    xlabel = "Frequency [cpd]",
    ylabel = "Φ [(m/s)²/cpd]",
    title = "Windspeed Spectra"
    )
plot!(freq_subsample, Φ_subsample, xaxis = :log, yaxis = :log, labels = "subsampled spectrum")

plot(freq, Φ,
    xaxis = :log,
    yaxis = :log,
    xlims = [freq_subsample[1], freq_subsample[end]],
    legend = :bottomright,
    labels = "standard spectrum",
    xlabel = "Frequency [cpd]",
    ylabel = "Φ [(m/s)²/cpd]",
    title = "Windspeed Spectra"
    )
plot!(freq_subsample, Φ_subsample, xaxis = :log, yaxis = :log, labels = "subsampled spectrum")


plot(freq, Φ,
    xaxis = :log,
    yaxis = :log,
    xlims = [1, 2],
    legend = :bottomright,
    labels = "standard spectrum",
    xlabel = "Frequency [cpd]",
    ylabel = "Φ [(m/s)²/cpd]",
    title = "Windspeed Spectra"
    )
plot!(freq_subsample, Φ_subsample, xaxis = :log, yaxis = :log, labels = "subsampled spectrum")
vline!([f_alias(freq_subsample[end], 1/(0.5))], labels = "semidiurnal alias freq.")


@info sum(Φ_subsample * freq_subsample[1]) / (sum(Φ) * freq[1])