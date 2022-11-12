export loadms

abstract type AbstractMeasurementSet{T} end

struct AntennaTable{T}
stationcode::Array{String,1}
mount::Array{String,1}
dishdiameter::Array{T,1}
end

struct CasacorejlMeasurementSet{T} <: AbstractMeasurementSet{T}
uvw::Array{T,2}
datavis::Vector{Array{Complex{T},2}}
antenna1::Array{Int,1}
antenna2::Array{Int,1}
times::Array{T,1}
scannumber::Array{Int,1}
weight::Vector{Array{T,1}}
weightspectrum::Array{T,3}
sigma::Vector{Array{T,1}}
sigmaspectrum::Array{T,3}
antennatable::AntennaTable{T}
end

struct CasacorePyMeasurementSet{T} <: AbstractMeasurementSet{T} end

function loadms(msname::String, stations::String)
    tab = CasacoreTable(msname, CasacoreTables.Old)

    # create MeasurementSet
    uvw = tab[:UVW][:, :] # uvw can be indexed with : (but other arrays need explicit limits)
    nrows = size(uvw)[2]
    datavis = tab[:DATA][1:nrows]
    antenna1 = tab[:ANTENNA1][:]
    antenna2 = tab[:ANTENNA2][:]
    times = tab[:TIME][:]
    scannumber = tab[:SCAN_NUMBER][:]
    weight = tab[:WEIGHT][1:nrows]
    weightspectrum = tab[:WEIGHT_SPECTRUM][:, :, :]
    sigma = tab[:SIGMA][1:nrows]
    sigmaspectrum = tab[:SIGMA_SPECTRUM][:, :, :]

    # create AntennaTable
    df = CSV.read(stations, DataFrame; delim=" ", ignorerepeated=true)
    antennatable = AntennaTable{Float64}(df.station,df.mount,df.dishdiameter_m)

    measurementset = CasacorejlMeasurementSet{Float64}(uvw,datavis,antenna1,antenna2,times,scannumber,weight,weightspectrum,sigma,sigmaspectrum,antennatable)

    return measurementset
end
