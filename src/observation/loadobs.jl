export loadobs

abstract type AbstractObservation{T} end

struct CjlObservation{T} <: AbstractObservation{T}
    uvw::Matrix{Float64}
    data::Array{Complex{Float32},4}
    antenna1::Vector{Int}
    antenna2::Vector{Int}
    times::Vector{Float64}
    exposure::Float64
    scanno::Vector{Int}
    weight::Vector{Vector{Float32}}
    weightspec::Array{Float32,3}
    sigma::Vector{Vector{Float32}}
    sigmaspec::Array{Float32,3}
    numchan::Int64
    chanfreqvec::Array{Float64,1}
    chanwidth::Float64
    stationinfo::DataFrame
    yamlconf::Dict
    rngcorrupt::AbstractRNG
    rngtrop::AbstractRNG
end

struct CpyObservation{T} <: AbstractObservation{T} end

#=struct StationInfo{T} <: Tables.AbstractColumns
    station::Array{String,1}
    dishdiameter::Array{T,1}
    sefd::Array{T,1}
    pwv::Array{T,1}
    groundpresssure::Array{T,1}
    groundtemperature::Array{T,1}
    coherencetime::Array{T,1}
    pointingrms::Array{T,1}
    beamfwhm230::Array{T,1}
    beammodel::Array{String,1}
    apefficiency::Array{T,1}
    g_pol1_re_mean::Array{T,1}
    g_pol1_re_std::Array{T,1}
    g_pol1_im_mean::Array{T,1}
    g_pol1_im_std::Array{T,1}
    g_pol2_re_mean::Array{T,1}
    g_pol2_re_std::Array{T,1}
    g_pol2_im_mean::Array{T,1}
    g_pol2_im_std::Array{T,1}
    d_pol1_re_mean::Array{T,1}
    d_pol1_re_std::Array{T,1}
    d_pol1_im_mean::Array{T,1}
    d_pol1_im_std::Array{T,1}
    d_pol2_re_mean::Array{T,1}
    d_pol2_re_std::Array{T,1}
    d_pol2_im_mean::Array{T,1}
    d_pol2_im_std::Array{T,1}
    feedangle::Array{T,1}
    mount::Array{String,1}
end=#

# implement Tables interface for CjlObservation
Tables.istable(::Type{<:CjlObservation}) = true
Tables.columnaccess(::Type{<:CjlObservation}) = true
Tables.columns(m::CjlObservation) = m

# implement Tables interface for StationInfo
#=Tables.istable(::Type{<:StationInfo}) = true
Tables.columnaccess(::Type{<:StationInfo}) = true
Tables.columns(m::StationInfo) = m=#

# define loadobs function
#function loadobs(msname::String, stationinfo::String, delim::String, ignorerepeated::Bool)
function loadobs(yamlconf::Dict, delim::String, ignorerepeated::Bool)
    """
    load data and metadata from ms and station table
    """
    tab = CCTable(yamlconf["msname"], CCTables.Old)

    # read values from ms
    uvw = tab[:UVW][:,:]
    data = tab[:DATA][:]
    antenna1 = tab[:ANTENNA1][:]
    antenna2 = tab[:ANTENNA2][:]
    times = tab[:TIME][:]
    exposure = tab[:EXPOSURE][1]
    scanno = tab[:SCAN_NUMBER][:]
    weight = tab[:WEIGHT][:]
    weightspec = tab[:WEIGHT_SPECTRUM][:,:,:]
    sigma = tab[:SIGMA][:]
    sigmaspec = tab[:SIGMA_SPECTRUM][:,:,:]

    spectab = tab.SPECTRAL_WINDOW
    numchan = spectab[:NUM_CHAN][1]
    chanfreqvec = spectab[:CHAN_FREQ][1]
    chanwidth = spectab[:CHAN_WIDTH][1][1] # get only the first element instead of the entire channel width vector

    # reshape the various arrays
    data3d = reduce((x,y) -> cat(x, y, dims=3), data)
    data3dres = reshape(data3d, 2, 2, size(data[1])[2], :) # get nchan as 3rd dim and all rows as 4th dim
    data3dresandperm = permutedims(data3dres, (2,1,3,4))

    #=flag3d = reduce((x,y) -> cat(x, y, dims=3), flag)
    flag3dres = reshape(flag3d, 2, 2, size(flag[1])[2], :) # get nchan as 3rd dim and all rows as 4th dim
    flag3dresandperm = permutedims(flag3dres, (2,1,3,4))=#

    #observation = CjlObservation{Float64}(uvw,data,antenna1,antenna2,times,scanno,weight,weightspec,sigma,sigmaspec)
    
    # construct ms dataframe
    #=measurementset = DataFrame()
    measurementset.uvw = uvw
    measurementset.data = data
    measurementset.antenna1 = antenna1
    measurementset.antenna2 = antenna2
    measurementset.times = times
    measurementset.scanno = scanno
    measurementset.weight = weight
    measurementset.weightspectrum = weightspec
    measurementset.sigma = sigma
    measurementset.sigmaspectrum = sigmaspec=#

    # read values from station info csv file
    stationinfo = CSV.read(yamlconf["stations"], DataFrame; delim=delim, ignorerepeated=ignorerepeated)

    # generate some quantities to be available for all corrupting functions and them to the observation composite type
    @info("Creating random number generator instance seeded with $(yamlconf["corruptseed"]) for all corruptions except troposphere...")
    rngcorrupt = Xoshiro(Int(yamlconf["corruptseed"]))

    @info("Creating random number generator instance seeded with $(yamlconf["troposphere"]["tropseed"]) for troposphere...")
    rngtrop = Xoshiro(Int(yamlconf["troposphere"]["tropseed"]))

    # construct CjlObservation object
    observation = CjlObservation{Float64}(uvw,data3dresandperm,antenna1,antenna2,times,exposure,scanno,weight,weightspec,sigma,sigmaspec,
					  numchan,chanfreqvec,chanwidth,stationinfo,yamlconf,rngcorrupt,rngtrop)

    return observation
end
