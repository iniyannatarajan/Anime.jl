export loadms, CjlObservation

abstract type AbstractObservation{T} end

"""
Container type for storing observation parameters and data.
"""
struct CjlObservation{T} <: AbstractObservation{T}
    """
    Name of the Measurement Set
    """
    msname::String
    """
    Complex visibility data
    """
    data::Array{Complex{Float32},4}
    """
    Flag array of the same dimensions as data
    """
    flag::Array{Bool,4}
    """
    Flag row array that holds FLAG_ROW from MS
    """
    flagrow::Array{Bool,1}
    """
    Antenna 1 in a baseline pair
    """
    antenna1::Vector{Int32}
    """
    Antenna 2 in a baseline pair
    """
    antenna2::Vector{Int32}
    """
    uvw coordinates
    """
    uvw::Matrix{Float64}
    """
    Timestamps of complex visibilities
    """
    times::Vector{Float64}
    """
    Integration time
    """
    exposure::Float64
    """
    Scan numbers from MS
    """
    scanno::Vector{Int32}
    #=weight::Vector{Vector{Float32}}
    weightspec::Array{Float32,3}
    sigma::Vector{Vector{Float32}}
    sigmaspec::Array{Float32,3}=#
    """
    Number of frequency channels
    """
    numchan::Int64
    """
    Channel frequencies (Hz)
    """
    chanfreqvec::Array{Float64,1}
    """
    Width of frequency channel
    """
    chanwidth::Float64
    """
    Direction of phase centre
    """
    phasedir::Array{Float64,2}
    """
    Antenna positions in x,y,z coordinates
    """
    pos::Array{Float64, 2}
    """
    Dataframe of station information input by user
    """
    stationinfo::DataFrame
    #yamlconf::Dict
    """
    Consider only water vapour in the troposphere
    """
    tropwetonly::Bool
    """
    Correlator efficiency
    """
    correff::Float64
    """
    Introduce attenuation by the troposphere
    """
    tropattenuate::Bool
    """
    Introduce noise due to troposphere
    """
    tropskynoise::Bool
    """
    Introduce mean delays due to troposphere
    """
    tropmeandelays::Bool
    """
    Introduce turbulence in troposphere
    """
    tropturbulence::Bool
    """
    Polarization frame
    """
    polframe::String
    """
    Polarization mode
    """
    polmode::String
    """
    Pointing interval (s)
    """
    ptginterval::Float64
    """
    Pointing error mode
    """
    ptgmode::String
    """
    Station gain mode
    """
    stationgainsmode::String
    """
    Input bandpass file
    """
    bandpassfile::String
    """
    RNG for all models except troposphere
    """
    rngcorrupt::AbstractRNG
    """
    RNG for tropospheric models
    """
    rngtrop::AbstractRNG
end

struct CpyObservation{T} <: AbstractObservation{T} end

# implement Tables interface for CjlObservation
#Tables.istable(::Type{<:CjlObservation}) = true
#Tables.columnaccess(::Type{<:CjlObservation}) = true
#Tables.columns(m::CjlObservation) = m

"""
    loadms(msname::String, stations::String, corruptseed::Int64, tropseed::Int64, tropwetonly::Bool, correff::Float64, tropattenuate::Bool,
    tropskynoise::Bool, tropmeandelays::Bool, tropturbulence::Bool, polframe::String, polmode::String, ptginterval::Float64, ptgmode::String,
    stationgainsmode::String, bandpassfile::String; delim::String=",", ignorerepeated::Bool=false)

Load data and metadata from MS, station and bandpass tables and return a CjlObservation object.
"""
function loadms(msname::String, stations::String, corruptseed::Int64, tropseed::Int64, tropwetonly::Bool, correff::Float64, tropattenuate::Bool,
    tropskynoise::Bool, tropmeandelays::Bool, tropturbulence::Bool, polframe::String, polmode::String, ptginterval::Float64, ptgmode::String,
    stationgainsmode::String, bandpassfile::String; delim::String=",", ignorerepeated::Bool=false)

    tab = CCTable(msname, CCTables.Old)

    # read values from ms
    data::Vector{Matrix{ComplexF32}} = tab[:DATA][:]
    flag::Vector{Matrix{Bool}} = tab[:FLAG][:]
    flagrow::Vector{Bool} = tab[:FLAG_ROW][:]
    antenna1::Vector{Int32} = tab[:ANTENNA1][:]
    antenna2::Vector{Int32} = tab[:ANTENNA2][:]
    uvw::Matrix{Float64} = tab[:UVW][:,:]
    times::Vector{Float64} = tab[:TIME][:]
    exposure::Float64 = tab[:EXPOSURE][1]
    scanno::Vector{Int32} = tab[:SCAN_NUMBER][:]
    #=weight::Vector{Vector{Float32}} = tab[:WEIGHT][:]
    weightspec::Array{Float32, 3} = tab[:WEIGHT_SPECTRUM][:,:,:]
    sigma::Vector{Vector{Float32}} = tab[:SIGMA][:]
    sigmaspec::Array{Float32, 3} = tab[:SIGMA_SPECTRUM][:,:,:]=#

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

    flag3d = reduce((x,y) -> cat(x, y, dims=3), flag)
    flag3dres = reshape(flag3d, 2, 2, size(flag[1])[2], :) # get nchan as 3rd dim and all rows as 4th dim
    flag3dresandperm = permutedims(flag3dres, (2,1,3,4))

    # read values from station info csv file
    stationinfo = CSV.read(stations, DataFrame; delim=delim, ignorerepeated=ignorerepeated)

    # parse strings to complex values for gjones terms
    stationinfo.g_pol1_loc = map(x->parse(ComplexF32,x), stationinfo.g_pol1_loc)
    stationinfo.g_pol2_loc = map(x->parse(ComplexF32,x), stationinfo.g_pol2_loc)
    stationinfo.g_pol1_scale = map(x->parse(ComplexF32,x), stationinfo.g_pol1_scale)
    stationinfo.g_pol2_scale = map(x->parse(ComplexF32,x), stationinfo.g_pol2_scale)

    # parse strings to complex values for djones terms
    stationinfo.d_pol1_loc = map(x->parse(ComplexF32,x), stationinfo.d_pol1_loc)
    stationinfo.d_pol2_loc = map(x->parse(ComplexF32,x), stationinfo.d_pol2_loc)
    stationinfo.d_pol1_scale = map(x->parse(ComplexF32,x), stationinfo.d_pol1_scale)
    stationinfo.d_pol2_scale = map(x->parse(ComplexF32,x), stationinfo.d_pol2_scale)

    # parse strings t
    stationinfo.pbmodel = map(strip, stationinfo.pbmodel)
    stationinfo.mount = map(strip, stationinfo.mount)

    # generate some quantities to be available for all corrupting functions and them to the observation composite type
    rngcorrupt = Xoshiro(corruptseed)
    rngtrop = Xoshiro(tropseed)

    # construct CjlObservation object    
    observation = CjlObservation{Float64}(msname,data3dresandperm,flag3dresandperm,flagrow,antenna1,antenna2,uvw,times,exposure,scanno,numchan,chanfreqvec,
    chanwidth,phasedir,pos,stationinfo,tropwetonly,correff,tropattenuate,tropskynoise,tropmeandelays,tropturbulence,polframe,polmode,
    ptginterval,ptgmode,stationgainsmode,bandpassfile,rngcorrupt,rngtrop)

    @info("Load data for processing 🙆")
    return observation
end