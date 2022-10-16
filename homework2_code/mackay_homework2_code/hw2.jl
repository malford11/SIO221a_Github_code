using Plots, StatsPlots
using NetCDF, Dates
using Statistics, StatsBase, LinearAlgebra, Distributions

DATA_PATH = joinpath(@__DIR__, "../data")
FIGURE_PATH = joinpath(@__DIR__, "figures")

pier_data_filepath(year::Int) = joinpath(DATA_PATH, "scripps_pier-$year.nc")

# ds_2008 = Dataset(pier_data_filepath(2008))
# ds_2018 = Dataset(pier_data_filepath(2018))

    # time = unix2datetime.(time) # converts based on 1970-01-01 00:00:00 UTC convention
# 1.
pier_data_2008 = pier_data_filepath(2008)
pier_data_2018 = pier_data_filepath(2018)

temperature_2008 = ncread(pier_data_2008, "temperature") # °C
temperature_2018 = ncread(pier_data_2018, "temperature")
pressure_2008 = ncread(pier_data_2008, "pressure") # dbar
pressure_2018 = ncread(pier_data_2018, "pressure")
time_2008 = unix2datetime.(ncread(pier_data_2008, "time"))
time_2018 = unix2datetime.(ncread(pier_data_2018, "time"))

plot(time_2008, temperature_2008, title = "2008 Pier Temperature", xlabel = "time [yyyy-mm-dd]", ylabel = "temperature [°C]")
savefig(joinpath(FIGURE_PATH, "temp_2008.png"))
plot(time_2018, temperature_2018, title = "2018 Pier Temperature", xlabel = "time [yyyy-mm-dd]", ylabel = "temperature [°C]")
savefig(joinpath(FIGURE_PATH, "temp_2018.png"))
plot(time_2008, pressure_2008, title = "2008 Pier Pressure", xlabel = "time [yyyy-mm-dd]", ylabel = "pressure [dbar]")
savefig(joinpath(FIGURE_PATH, "pres_2008.png"))
plot(time_2018, pressure_2018, title = "2018 Pier Pressure", xlabel = "time [yyyy-mm-dd]", ylabel = "pressure [dbar]")
savefig(joinpath(FIGURE_PATH, "pres_2018.png"))

# 2.
@info string("2008 Temperature: ", mean(temperature_2008), " ± ", std(temperature_2008))
@info string("2018 Temperature: ", mean(temperature_2018), " ± ", std(temperature_2018))
@info string("2008 Pressure: ", mean(pressure_2008), " ± ", std(pressure_2008))
@info string("2018 Pressure: ", mean(pressure_2018), " ± ", std(pressure_2018))

# subsample
daily_temp_2008 = temperature_2008[1:180:end] # sampled every 8 min
daily_temp_2018 = temperature_2018[1:360:end] # sampled every 4 min
daily_pres_2008 = pressure_2008[1:180:end]
daily_pres_2018 = pressure_2018[1:360:end]

@info string("2008 Daily Temperature: ", mean(daily_temp_2008), " ± ", std(daily_temp_2008))
@info string("2018 Daily Temperature: ", mean(daily_temp_2018), " ± ", std(daily_temp_2018))
@info string("2008 Daily Pressure: ", mean(daily_pres_2008), " ± ", std(daily_pres_2008))
@info string("2018 Daily Pressure: ", mean(daily_pres_2018), " ± ", std(daily_pres_2018))

# 3.
std_err(data) = var(data) * sqrt(2 / (length(data) - 1))
@info string("2008 Temperature: ", var(temperature_2008), " ± ", std_err(temperature_2008))
@info string("2018 Temperature: ", var(temperature_2018), " ± ", std_err(temperature_2018))
@info string("2008 Pressure: ", var(pressure_2008), " ± ", std_err(pressure_2008))
@info string("2018 Pressure: ", var(pressure_2018), " ± ", std_err(pressure_2018))

@info string("2008 Daily Temperature: ", var(daily_temp_2008), " ± ", std_err(daily_temp_2008))
@info string("2018 Daily Temperature: ", var(daily_temp_2018), " ± ", std_err(daily_temp_2018))
@info string("2008 Daily Pressure: ", var(daily_pres_2008), " ± ", std_err(daily_pres_2008))
@info string("2018 Daily Pressure: ", var(daily_pres_2018), " ± ", std_err(daily_pres_2018))

# 4.
# 2008
density(temperature_2008, legend = false)
savefig(joinpath(FIGURE_PATH, "epdf_2008.png"))
std2 = mean(temperature_2008) + 2 * std(temperature_2008)
@info sum(temperature_2008 .>= std2) / length(temperature_2008)

# 2018
density(temperature_2018, legend = false)
savefig(joinpath(FIGURE_PATH, "epdf_2018.png"))
std2 = mean(temperature_2018) + 2 * std(temperature_2018)
@info sum(temperature_2018 .>= std2) / length(temperature_2018)

# 5.
# temperature
density(temperature_2008, label = "2008 epdf", linecolor = :blue, linewidth = 2)
density!(temperature_2018, label = "2018 epdf", linecolor = :red, linewidth = 2)
plot!(
    Normal(mean(temperature_2008),
    std(temperature_2008)),
    label = "2008 Gaussian",
    linestyle = :dash,
    linecolor = :blue,
    linewidth = 2
)
plot!(
    Normal(mean(temperature_2018),
    std(temperature_2018)),
    label = "2018 Gaussian",
    linestyle = :dash,
    linecolor = :red,
    linewidth = 2
)
bounds_2008 = [-1 1; 1 1] \ [2 * sqrt(3) * std(temperature_2008); 2 * mean(temperature_2008)]
bounds_2018 = [-1 1; 1 1] \ [2 * sqrt(3) * std(temperature_2018); 2 * mean(temperature_2018)]
plot!(
    Uniform(bounds_2008[1], bounds_2008[2]),
    label = "2008 Uniform",
)
plot!(
    Uniform(bounds_2018[1], bounds_2018[2]),
    label = "2018 Uniform",
)
title!("Temperature PDFs")
xlabel!("Temperature [°C]")
savefig(joinpath(FIGURE_PATH, "epdf_temp.png"))

# pressure
density(pressure_2008, label = "2008 epdf", linecolor = :blue, linewidth = 2)
density!(pressure_2018, label = "2018 epdf", linecolor = :red, linewidth = 2)
plot!(
    Normal(mean(pressure_2008),
    std(pressure_2008)),
    label = "2008 Gaussian",
    linestyle = :dash,
    linecolor = :blue,
    linewidth = 2
)
plot!(
    Normal(mean(pressure_2018),
    std(pressure_2018)),
    label = "2018 Gaussian",
    linestyle = :dash,
    linecolor = :red,
    linewidth = 2
)
# for a uniformly distributed random variable
# μ = (a + b) / 2, σ² = 1/12(b-a)²
bounds_2008 = [-1 1; 1 1] \ [2 * sqrt(3) * std(pressure_2008); 2 * mean(pressure_2008)]
bounds_2018 = [-1 1; 1 1] \ [2 * sqrt(3) * std(pressure_2018); 2 * mean(pressure_2018)]
plot!(
    Uniform(bounds_2008[1], bounds_2008[2]),
    label = "2008 Uniform",
)
plot!(
    Uniform(bounds_2018[1], bounds_2018[2]),
    label = "2018 Uniform",
)
title!("Pressure PDFs")
xlabel!("Pressure [dbars]")
savefig(joinpath(FIGURE_PATH, "epdf_pres.png"))

# 6.
k = [1, 45, 90, 180, 360, 720, 1440, 2880, 5760]
x = zeros(5)
N = zeros(5)
for i in 1:5
    data = temperature_2018[1:k[i]:end]
    N[i] = length(data)
    x[i] = var(data)
end
plot(
    N,
    x,
    xlabel = "Inverse Sample Size (1/N)",
    ylabel = "Sample Variance (°C)²",
    title = "Temperature N-Sample Variance "
)
savefig(joinpath(FIGURE_PATH, "nsample.png"))
