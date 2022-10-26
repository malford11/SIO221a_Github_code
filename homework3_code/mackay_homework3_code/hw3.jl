using CSVFiles, DataFrames
using Dates, Statistics, LinearAlgebra
using Plots

DATA_PATH = joinpath(@__DIR__, "data/")
FIG_PATH = joinpath(@__DIR__, "figures/")

data2015, data2016 = map((2015, 2016)) do yr
    filename = joinpath(DATA_PATH, "46047h$yr.txt.gz")
    data = DataFrame(
            load(
                File(format"CSV", filename);
                spacedelim = true,
                header_exists = false,
                skiplines_begin = 2,
                colnames = ["YY", "MM", "DD", "hh", "mm", "WDIR", "WSPD", "GST",
                            "WVHT", "DPD", "APD", "MWD", "PRES", "ATMP", "WTMP", "DEWP", "VIS", "TIDE"] # yr mo dy hr mn degT m/s m/s m sec sec degT hPa degC degC degC mi ft
            )
        )
    # Create single DateTime
    transform!(data, [:YY, :MM, :DD, :hh, :mm] => ByRow((YY, MM, DD, hh, mm) -> DateTime(YY, MM, DD, hh, mm)) => :Time)

    # Replace missing values with NaN
    replace!(data.WSPD, 99 => NaN)
    replace!(data.WVHT, 99 => NaN)
    replace!(data.WTMP, 999 => NaN)
    replace!(data.ATMP, 999 => NaN)
    return data
end
data = vcat(data2015, data2016)

# 1. Plot the time series of wind speed, wave height, water temperature, and air temperature from 46047
plot(data.Time, data.WSPD, ylabel = "Wave Speed [m/s]", xlabel = "Time [YYYY-MM-DD]", title = "2015-16 Wave Speed")
savefig(joinpath(FIG_PATH, "wspd.png"))
plot(data.Time, data.WVHT, ylabel = "Wave Height [m]", xlabel = "Time [YYYY-MM-DD]", title = "2015-16 Wave Height")
savefig(joinpath(FIG_PATH, "wvht.png"))
plot(data.Time, data.WTMP, ylabel = "Sea Surface Temp [°C]", xlabel = "Time [YYYY-MM-DD]", title = "2015-16 SST")
savefig(joinpath(FIG_PATH, "wtmp.png"))
plot(data.Time, data.ATMP, ylabel = "Air Temp [°C]", xlabel = "Time [YYYY-MM-DD]", title = "2015-16 Air Temperature")
savefig(joinpath(FIG_PATH, "atmp.png"))

# 2. Average the data to produce monthly means for 2015 and 2016. Plot the means for each month and standard error of the mean.
function monthly_mean(data::DataFrame, col::Symbol, mm::Int, yy::Int)
    mean(data[!,col][month.(data.Time) .== mm .&& year.(data.Time) .== yy .&& .!isnan.(data[!,col])])
end

function monthly_err(data::DataFrame, col::Symbol, mm::Int, yy::Int)
    data = data[!,col][month.(data.Time) .== mm .&& year.(data.Time) .== yy .&& .!isnan.(data[!,col])]
    subdata = data[1:7*24:end] #subsample every 7 days
    return var(subdata) / sqrt(length(subdata))
end

function monthly_summary(data::DataFrame, col::Symbol)
    x = zeros(12, 2)
    y = similar(x)
    t = similar(x, DateTime)
    for (j, yy) in enumerate((2015, 2016))
        for (i, mm) in enumerate(1:12)
            x[i, j] = monthly_mean(data, col, mm, yy)
            y[i, j] = monthly_err(data, col, mm, yy)
            t[i, j] = DateTime(yy, mm)
        end
    end
    return vec(t), vec(x), vec(y)
end

x, y, e = monthly_summary(data, :WSPD)
plot(x, y, yerr = e, ylabel = "Wave Speed [m/s]", xlabel = "Time [YYYY-MM-DD]", title = "Monthly Mean Wave Speed", legend = false)
savefig(joinpath(FIG_PATH, "mean_wspd.png"))
x, y, e = monthly_summary(data, :WVHT)
plot(x, y, yerr = e, ylabel = "Wave Height [m]", xlabel = "Time [YYYY-MM-DD]", title = "Monthly Mean Wave Height", legend = false)
savefig(joinpath(FIG_PATH, "mean_wvht.png"))
x, y, e = monthly_summary(data, :WTMP)
savefig(joinpath(FIG_PATH, "mean_wtmp.png"))
plot(x, y, yerr = e, ylabel = "Sea Surface Temp [°C]", xlabel = "Time [YYYY-MM-DD]", title = "Monthly Mean SST", legend = false)
x, y, e = monthly_summary(data, :ATMP)
plot(x, y, yerr = e, ylabel = "Air Temp [°C]", xlabel = "Time [YYYY-MM-DD]", title = "Monthly Mean SST", legend = false)
savefig(joinpath(FIG_PATH, "mean_atmp.png"))

# 3. Treating the two years separately, least-squares fit a mean and an annual cycle to the four data records.
# 4. Augment your annual cycle least-squares fit with a semi-annual cycle.
#=
Note that we can forgo using the annual fitting function and just use the semi-annual fit because
of the orthogonality of the basis functions - the coefficients for the annual cycle will be unaffected
by the additional semiannual terms.
=#
function fit_annual_cycle(data::DataFrame, col::Symbol)
    filter = .!isnan.(data[!, col])
    ω = 2 * π / (365 * 24) # hr/yr
    t = 1:length(data.Time[filter])
    A = [ones(length(t)) sin.(ω .* t) cos.(ω .* t)]
    b = data[!, col][filter]
    return A \ b
end

function fit_semiannual_cycle(data::DataFrame, col::Symbol)
    filter = .!isnan.(data[!, col])
    ω1 = 2 * π / (365 * 24) # annual freq
    ω2 = 2 * π / (365 * 12) # semi-annual freq
    t = 1:length(data.Time[filter])
    A = [ones(length(t)) sin.(ω1 .* t) cos.(ω1 .* t) sin.(ω2 .* t) cos.(ω2 .* t)]
    b = data[!, col][filter]
    return A \ b
end

cycle_mean(coefs::Vector) = coefs[1]
annual_cycle_amplitude(coefs::Vector) = norm(coefs[2:3])
semiannual_cycle_amplitude(coefs::Vector) = norm(coefs[4:5])

# 5. What is the squared misfit of your least-squares fits?
function chi2_misfit(data, col, coefs)
    ω1 = 2 * π / (365 * 24) # annual freq
    ω2 = 2 * π / (365 * 12) # semi-annual freq

    filter = .!isnan.(data[!, col])
    t = 1:length(data.Time[filter])
    A = [ones(length(t)) sin.(ω1 .* t) cos.(ω1 .* t) sin.(ω2 .* t) cos.(ω2 .* t)]
    σ2 = var(data[!, col][filter][1:7*24:end])
    semiannual_chi2 = sum((data[!, col][filter] .- A * coefs).^2 ./ σ2)
    annual_chi2 = sum((data[!, col][filter] .- A[:, 1:3] * coefs[1:3]).^2 ./ σ2)
    return annual_chi2, semiannual_chi2
end


function least_squares_summary(data)
    colnames = (:C1, :C2, :C3, :C4, :C5, :MEAN, :ANNUAL_AMP, :SEMIANNUAL_AMP, :ANNUAL_CHI2, :SEMIANNUAL_CHI2)
    coltypes = map(x -> Float64, colnames)
    df = DataFrame(NamedTuple{colnames}(type[] for type in coltypes))
    df.VARNAME = Symbol[]
    for var in (:WSPD, :WVHT, :WTMP, :ATMP)
        coefs = fit_semiannual_cycle(data, var)
        push!(df, (Tuple(coefs)..., cycle_mean(coefs), annual_cycle_amplitude(coefs), semiannual_cycle_amplitude(coefs), chi2_misfit(data, var, coefs)..., var))
    end
    return df
end

least_squares_summary(data2015)
least_squares_summary(data2016)
