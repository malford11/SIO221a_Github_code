#=
A script to auto-download SCCOOS SIO Pier data files
=#

using RemoteFiles, NetCDF

# where the data is saved
DATA_PATH = joinpath(@__DIR__, "data")

"""
    pier_url(year)

Returns the download URL for Scripps Pier data for a given year.
"""
pier_url(year::Int) = "https://thredds.sccoos.org/thredds/fileServer/autoss/scripps_pier-$year.nc"

# Create a set of all files to download
years = Tuple(2005:2021)
files = map(years) do year
    Symbol("pier_data_$year"), RemoteFile(pier_url(year), dir = DATA_PATH)
end
pier_data = RemoteFileSet("SCCOOS Scripps pier sensor data"; files...)

download(pier_data)