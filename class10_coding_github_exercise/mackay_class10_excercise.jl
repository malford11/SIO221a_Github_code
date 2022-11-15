using MAT, Test
using Statistics, FastTransforms
using Plots

include("../SpectralAnalysis.jl") # include SpectralAnalysis module - functions live here

# read matlab file - not worth creating a function since the internal dict
# storage may change from project to project making it difficult to standardize
file = joinpath(@__DIR__, "class10_record.mat")
vars = matread(file)
data = vars["class10_record"]["data"]' # u vel in m/s from 300 m, station PAPA.
time = vars["class10_record"]["time"]' # Time is datenum [days]; sample interval is 30 min. 
# A serial date number represents the whole and fractional number of days from a fixed, preset date 

T = (time[end] - time[1]) / 10 # create 10 segments
Φ, f, (σ_low, σ_high) = SpectralAnalysis.spectrum(data, time, T)

plot(f, Φ,
    xlabel = "frequency [cpd]",
    ylabel = "Spectral Density [(m/s)² / cpd]",
    title = "Zonal Velocity Spectrum at 300m\nStation PAPA",
    ribbon = (Φ .* σ_low, Φ .* σ_high),
    legend = false,
    xaxis = :log,
    yaxis = :log
)

# TESTS
@testset "Demean & Detrend" begin
    # test demean of data gives data with zero mean
    @test isapprox(0.0, mean(SpectralAnalysis.demean(data)); atol = sqrt(eps(Float64)))
    # trend of detrend gives zero y-intercept & zero slope
    @test isapprox([0.0; 0.0], SpectralAnalysis.trend(SpectralAnalysis.detrend(data, time), time); atol = sqrt(eps(Float64)))
end

@testset "Segmenting" begin
    a = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    @test size(SpectralAnalysis.segment(a, 2)) == (2, 9)
    @test size(SpectralAnalysis.segment(a, 4)) == (4, 4)
    @test size(SpectralAnalysis.segment(a, 6)) == (6, 2)
end

@testset "Spectrum" begin
    # check spectrum function gives expected result for single segment
    # a better test is probably to check the spectrum for something we can compare to analytically, like a sine wave
    X1 = fft(SpectralAnalysis.detrend(data, time)) ./ length(data)
    Φ1 = (norm.(X1).^2 / (1 / (time[end] - time[1])))[1:Int(length(data)/2)] # expected
    dt = time[2] - time[1]
    T = length(data) * dt
    Φ2, _, _ = SpectralAnalysis.spectrum(data, time, T) # function result
    @test isapprox(Φ1, Φ2)
end
