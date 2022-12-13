export loadobs

abstract type AbstractObservation{T} end

# TODO the entire data structure to hold MS data needs to be rewritten
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
    phasedir::Array{Float64,2}
    pos::Array{Float64, 2}
    stationinfo::DataFrame
    yamlconf::Dict
    rngcorrupt::AbstractRNG
    rngtrop::AbstractRNG
end

struct CpyObservation{T} <: AbstractObservation{T} end

# implement Tables interface for CjlObservation
Tables.istable(::Type{<:CjlObservation}) = true
Tables.columnaccess(::Type{<:CjlObservation}) = true
Tables.columns(m::CjlObservation) = m

# define loadobs function
function loadobs(yamlconf::Dict, delim::String, ignorerepeated::Bool)
    """
    load data and metadata from ms and station table
    """
    tab = CCTable(yamlconf["msname"], CCTables.Old)

    # read values from ms
    uvw::Matrix{Float64} = tab[:UVW][:,:]
    data::Vector{Matrix{ComplexF32}} = tab[:DATA][:]
    antenna1::Vector{Int32} = tab[:ANTENNA1][:]
    antenna2::Vector{Int32} = tab[:ANTENNA2][:]
    times::Vector{Float64} = tab[:TIME][:]
    exposure::Float64 = tab[:EXPOSURE][1]
    scanno::Vector{Int32} = tab[:SCAN_NUMBER][:]
    weight::Vector{Vector{Float32}} = tab[:WEIGHT][:]
    weightspec::Array{Float32, 3} = tab[:WEIGHT_SPECTRUM][:,:,:]
    sigma::Vector{Vector{Float32}} = tab[:SIGMA][:]
    sigmaspec::Array{Float32, 3} = tab[:SIGMA_SPECTRUM][:,:,:]

    spectab = tab.SPECTRAL_WINDOW
    numchan::Int32 = spectab[:NUM_CHAN][1]
    chanfreqvec::Vector{Float64} = spectab[:CHAN_FREQ][1]
    chanwidth::Float64 = spectab[:CHAN_WIDTH][1][1] # get only the first element instead of the entire channel width vector

    fieldtab = tab.FIELD
    phasedir::Matrix{Float64} = fieldtab[:PHASE_DIR][:][1] # 2 x 1 Matrix of ra and dec

    anttab = tab.ANTENNA
    pos::Matrix{Float64} = anttab[:POSITION][:,:]

    # reshape the various arrays
    data3d = reduce((x,y) -> cat(x, y, dims=3), data)
    data3dres = reshape(data3d, 2, 2, size(data[1])[2], :) # get nchan as 3rd dim and all rows as 4th dim
    data3dresandperm = permutedims(data3dres, (2,1,3,4))

    #=flag3d = reduce((x,y) -> cat(x, y, dims=3), flag)
    flag3dres = reshape(flag3d, 2, 2, size(flag[1])[2], :) # get nchan as 3rd dim and all rows as 4th dim
    flag3dresandperm = permutedims(flag3dres, (2,1,3,4))=#

    # read values from station info csv file
    stationinfo = CSV.read(yamlconf["stations"], DataFrame; delim=delim, ignorerepeated=ignorerepeated)

    # parse strings to complex values for gjones terms
    stationinfo.g_pol1_loc = map(x->parse(ComplexF32,x), stationinfo.g_pol1_loc)
    stationinfo.g_pol2_loc = map(x->parse(ComplexF32,x), stationinfo.g_pol2_loc)

    # parse strings to complex values for djones terms
    stationinfo.d_pol1_loc = map(x->parse(ComplexF32,x), stationinfo.d_pol1_loc)
    stationinfo.d_pol2_loc = map(x->parse(ComplexF32,x), stationinfo.d_pol2_loc)

    # parse strings t
    stationinfo.pbmodel = map(x->strip(x), stationinfo.pbmodel)
    stationinfo.mount = map(x->strip(x), stationinfo.mount)

    # generate some quantities to be available for all corrupting functions and them to the observation composite type
    rngcorrupt = Xoshiro(Int(yamlconf["corruptseed"]))
    rngtrop = Xoshiro(Int(yamlconf["troposphere"]["tropseed"]))

    # construct CjlObservation object
    observation = CjlObservation{Float64}(uvw,data3dresandperm,antenna1,antenna2,times,exposure,scanno,weight,weightspec,sigma,sigmaspec,
					  numchan,chanfreqvec,chanwidth,phasedir,pos,stationinfo,yamlconf,rngcorrupt,rngtrop)

    @info("Load observation and metadata into memory for processing... 🙆")
    return observation
end
