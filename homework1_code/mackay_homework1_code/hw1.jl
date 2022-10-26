using NetCDF
using Statistics, StatsBase, LinearAlgebra
using Plots, Dates

DATA_PATH = joinpath(@__DIR__, "../data")
FIGURE_PATH = joinpath(@__DIR__, "figures")

pier_data_filepath(year::Int) = joinpath(DATA_PATH, "scripps_pier-$year.nc")

# Problem 1 - 2021 Pier Analysis
# read pier data
pier_data_2021 = pier_data_filepath(2021)
temperature = ncread(pier_data_2021, "temperature") # °C
time = ncread(pier_data_2021, "time") # seconds since 1970-01-01 00:00:00 UTC
time = unix2datetime.(time) # converts based on 1970-01-01 00:00:00 UTC convention

plot(time,
    temperature; 
    plot_title = "Seawater Temperature at Scripps Pier\n2021",
    xlabel = "Time yyyy-mm-dd",
    ylabel = "Temperature °C",
    legend = false,
    )
savefig(joinpath(FIGURE_PATH, "temp_2021.png"))

# (μ = 17.273985f0, σ = 2.4253728f0)
@info (μ = mean(temperature), σ = std(temperature))

h = fit(Histogram, temperature)
h = normalize(h, mode = :pdf)
plot(h;
    plot_title = "Emperical pdf of Seawater Temperature\n2021",
    xlabel = "Temperate °C",
    legend = false,
    )
savefig(joinpath(FIGURE_PATH, "epdf_2021.png"))


# Problem 2 - 2005-2021 Pier Analysis
years = Tuple(2005:2021)
temperature = map(years) do year
    temperature = ncread(pier_data_filepath(year), "temperature") # °C
end
temperature = vcat(temperature...) # combines data for all years into single vector

time = map(years) do year
    time = ncread(pier_data_filepath(year), "time") # seconds since 1970-01-01 00:00:00 UTC
    time = unix2datetime.(time) # converts based on 1970-01-01 00:00:00 UTC convention
end
time = vcat(time...)

filter = temperature .< 40 # filter out temperatures greater than 40°C
time, temperature = time[filter], temperature[filter]

plot(time,
    temperature;
    plot_title = "Seawater Temperature at Scripps Pier\n2005-2021",
    xlabel = "Time yyyy-mm-dd",
    ylabel = "Temperature °C",
    legend = false,
    )
savefig(joinpath(FIGURE_PATH, "temp_2005_2021.png"))

# (μ = 17.799109f0, σ = 2.7387943f0)
@info (μ = mean(temperature), σ = std(temperature))

h = fit(Histogram, temperature)
h = normalize(h, mode = :pdf)
plot(h;
    plot_title = "Emperical pdf of Seawater Temperature\n2005-2021",
    xlabel = "Temperate °C",
    legend = false,
    )
savefig(joinpath(FIGURE_PATH, "ecdf_2005_2021.png"))
