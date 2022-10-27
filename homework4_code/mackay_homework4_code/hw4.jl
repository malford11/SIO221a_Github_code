using NetCDF, Dates
using Statistics, LinearAlgebra
using FastTransforms
using Plots

DATA_PATH = joinpath(@__DIR__, "data")
FIGURE_PATH = joinpath(@__DIR__, "figures")

pier_data_2021 = joinpath(DATA_PATH, "scripps_pier-2021.nc")
pressure = ncread(pier_data_2021, "pressure") # dbar
time = ncread(pier_data_2021, "time") # seconds since 1970-01-01T00:00:00

plot(unix2datetime.(time), pressure,
    title = "2021 Pier Pressure",
    xlabel = "Time [YYYY-MM-DD]",
    ylabel = "Pressure [dbar]",
    legend = false
)
savefig(joinpath(FIGURE_PATH, "pressure_2021.png"))

# fill gaps
startline=70521
endline=88190
xx=startline:endline
N=(time[endline] - time[startline]) / 240 # sample every 4 min = 240 s
time2 = time[startline] .+ 240 * (0:N)
pressure2 = similar(time2)
# gaps at 2420, 4952, 7496, 15042
pressure2[1:2420] = pressure[xx[1:2420]]
pressure2[2421] = mean(pressure[xx[2420:2421]])
pressure2[2422:4953] = pressure[xx[2421:4952]]
pressure2[4954] = mean(pressure[xx[4952:4953]])
pressure2[4955:7498] = pressure[xx[4953:7496]]
pressure2[7499] = mean(pressure[xx[7496:7497]])
pressure2[7500:15045] = pressure[xx[7497:15042]]
pressure2[15046] = mean(pressure[xx[15042:15043]])
pressure2[15047:length(xx)+4] = pressure[xx[15043:end]]

plot(unix2datetime.(time[2:end]), diff(time), ylims = (0, 1000),
    title = "Pier Data Measurement Time Increments",
    xlabel = "Time [YYYY-MM-DD]",
    ylabel = "Measurement Time Difference [s]",
    legend = false,
)
vline!(unix2datetime.([time[startline], time[endline]]))
savefig(joinpath(FIGURE_PATH, "timediff.png"))

# Least Squares Fit
## Periods [s]
O1 = 1 / (25.82 * 60 * 60) # principal lunar diurnal
K1 = 1 / (23.93 * 60 * 60) # luni-solar diurnal
M2 = 1 / (12.42 * 60 * 60) # principal lunar
A = [ones(length(time2)) sin.(2π * O1 .* time2) cos.(2π * O1 .* time2) sin.(2π * K1 .* time2) cos.(2π * K1 .* time2) sin.(2π * M2 .* time2) cos.(2π * M2 .* time2)]
coef = A \ pressure2
plot(unix2datetime.(time2), A * coef,
    title = "2021 Least-squares fit",
    xlabel = "Time [YYYY-MM-DD]",
    ylabel = "Pressure [dbar]",
    labels = "Least-squares"
)
plot!(unix2datetime.(time2), pressure2, labels = "Filled Data")
savefig(joinpath(FIGURE_PATH, "least_squares.png"))

# find amplitudes & mean:
@info "______LEAST SQUARES______"
@info "Mean (dbar): $(coef[1])"
@info "Amplitude O1 (dbar): $(norm(coef[2:3]))"
@info "Amplitude K1 (dbar): $(norm(coef[4:5]))"
@info "Amplitude M2 (dbar): $(norm(coef[6:7]))"

# Fourier Transform
N = length(pressure2)
data = pressure2
X_k = fftshift(fft(data)) ./ N

# Determine the mean
@info "______FOURIER TRANSFORM______"
@info "Mean: $(norm(fftshift(X_k)[1]))"

# Remove mean from further analysis
data = pressure2 .- mean(pressure2)
X_k = fftshift(fft(data)) ./ N

T = time2[end] - time2[1]
Δf = 1 / T
fN = (0.5 * N * Δf)
f = (-fN:Δf:fN)[1:end-1]

amp = 2 .* norm.(X_k)
amp_max = amp[amp .> 0.17]
f_max = f[amp .> 0.17]
for i in 1:3
    @info "Freq (1/s): $(abs(f_max[i])), Amplitude (dbar): $(amp_max[i]), Period (hr): $(1 / abs(f_max[i]) / 3600)"
end

plot(f, real.(X_k), xlims = (-0.00005, 0.00005), title = "Fourier Coefficients", xlabel = "Frequency [1/s]", labels = "Re(Xₖ)")
plot!(f, imag.(X_k), xlims = (-0.00005, 0.00005), xlabel = "Frequency [1/s]", labels = "Im(Xₖ)")
savefig(joinpath(FIGURE_PATH, "fourier_coef_components.png"))

plot(f, amp, xlims = (0, 0.00005), ylims = (0, 0.5),
     title = "Fourier Coefficient Amplitudes",
     xlabel = "Frequency [1/s]",
     ylabel = "|Xₖ| (dbar)",
     labels = "|Xₖ|")
vline!([O1], lw = 1, labels = "O1")
vline!([K1], lw = 1, labels = "K1")
vline!([M2], lw = 1, labels = "M2")
savefig(joinpath(FIGURE_PATH, "fourier_amplitudes.png"))

Φ = norm.(X_k).^2 ./ Δf
plot(f[Int((N/2)):end], Φ[Int((N/2)):end], title = "Spectral Energy Φ", xlabel = "Frequency [1/s]", legend = false)
savefig(joinpath(FIGURE_PATH, "spectral_energy.png"))
