export readms, readobsconfig, readstationinfo, readbandpassinfo, readalistv5, readalistv6, MeasurementSet

"""
    $(TYPEDEF)

Composite type for storing data from a Measurement Set.

# Fields
$(FIELDS)
"""
struct MeasurementSet{T}
    """
    Measurement Set name
    """
    name::String
    """
    Complex visibility data (MS `DATA` column)
    """
    data::Array{Complex{Float32},4}
    """
    Flag array of the same dimensions as data (MS `FLAG` column)
    """
    flag::Array{Bool,4}
    """
    Flag row Boolean vector (MS ```FLAG_ROW``` column)
    """
    flagrow::Array{Bool,1}
    """
    Antenna 1 in a baseline pair (MS `ANTENNA1` column)
    """
    antenna1::Vector{Int32}
    """
    Antenna 2 in a baseline pair (MS `ANTENNA2` column)
    """
    antenna2::Vector{Int32}
    """
    uvw coordinates (MS `UVW` column)
    """
    uvw::Matrix{Float64}
    """
    Timestamps (MS `TIME` column)
    """
    times::Vector{Float64}
    """
    Integration time
    """
    exposure::Float64
    """
    Scan numbers (```SCAN_NUMBER``` column in MS)
    """
    scanno::Vector{Int32}
    """
    Number of frequency channels
    """
    numchan::Int64
    """
    Channel frequencies (Hz) (MS ```CHAN_FREQ``` column)
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
end

"""
    readms(msname::String)

Read data from Measurement Set.
"""
function readms(msname::String)
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
    if :WEIGHT_SPECTRUM in keys(tab)
        weight_spectrum::Vector{Matrix{Float32}} = tab[:WEIGHT_SPECTRUM][:,:,:]
    else
        weight_spectrum = fill(fill(missing, (2, 2)), size(data))
    end

    sigma::Vector{Vector{Float32}} = tab[:SIGMA][:]
    if :SIGMA_SPECTRUM in keys(tab)
        sigma_spectrum::Vector{Matrix{Float32}} = tab[:SIGMA_SPECTRUM][:,:,:]
    else
        sigma_spectrum = fill(fill(missing, (2, 2)), size(data))
    end=#

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

    # populate MeasurementSet    
    measurementset = MeasurementSet{Float64}(msname, data3dresandperm, flag3dresandperm, flagrow, antenna1, 
    antenna2, uvw, times, exposure, scanno, numchan, chanfreqvec, chanwidth, phasedir, pos)

    @info("Load data from MS 🙆")
    return measurementset
end

"""
    readobsconfig(configfile::String)

Read YAML configuration file with observation settings in a Dictionary.
"""
function readobsconfig(configfile::String)
    obsconfig = YAML.load_file(configfile, dicttype=Dict{String,Any}) # load YAML config file

    return obsconfig
end

"""
    readstationinfo(stationfile::String; delim::String=",", ignorerepeated::Bool=false)

Read additional information about stations from CSV file into a DataFrame.
"""
function readstationinfo(stationfile::String; delim::String=",", ignorerepeated::Bool=false)
    # read values from station info csv file
    stationinfo = CSV.read(stationfile, DataFrame; delim=delim, ignorerepeated=ignorerepeated)

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

    # strip strings of superflous whitespaces
    stationinfo.pbmodel = map(strip, stationinfo.pbmodel)
    stationinfo.mount = map(strip, stationinfo.mount)
    
    return stationinfo
end

"""
    readbandpassinfo(bandpassfile::String; delim::String=",", ignorerepeated::Bool=false)

Read bandpass information from CSV file into a DataFrame.
"""
function readbandpassinfo(bandpassfile::String; delim::String=",", ignorerepeated::Bool=false)
    bandpassinfo = CSV.read(bandpassfile, DataFrame; delim=",", ignorerepeated=false)

    return bandpassinfo
end

"""
    readalistv5(alistfile::String)

Read in an alist v5 file generated by HOPS. 
"""
function readalistv5(alistfile::String)

    # set header names
    columns = "version,root_id,two,extent_no,duration,length,offset,expt_no,scan_id,procdate,year,timetag,scan_offset,source,baseline,quality,freq_code,polarization,lags,amp,snr,resid_phas,phase_snr,datatype,sbdelay,mbdelay,ambiguity,delay_rate,ref_elev,rem_elev,ref_az,rem_az,u,v,esdesp,epoch,ref_freq,total_phas,total_rate,total_mbdelay,total_sbresid,srch_cotime,noloss_cotime"
    header = [string(sub) for sub in split(columns, ",")]

    df = CSV.read(alistfile, DataFrame; header=header, comment="*", delim=" ", ignorerepeated=true)

    return df
end

"""
    readalistv6(alistfile::String)

Read in an alist v6 file generated by HOPS. 
"""
function readalistv6(alistfile::String)

    # set header names
    columns = "version,root_id,two,extent_no,duration,length,offset,expt_no,scan_id,procdate,year,timetag,scan_offset,source,baseline,quality,freq_code,polarization,lags,amp,snr,resid_phas,phase_snr,datatype,sbdelay,mbdelay,ambiguity,delay_rate,ref_elev,rem_elev,ref_az,rem_az,u,v,esdesp,epoch,ref_freq,total_phas,total_rate,total_mbdelay,total_sbresid,srch_cotime,noloss_cotime,ra_hrs,dec_deg,resid_delay"
    header = [string(sub) for sub in split(columns, ",")]

    df = CSV.read(alistfile, DataFrame; header=header, comment="*", delim=" ", ignorerepeated=true)

    return df
end

#=
    # generate some quantities to be available for all corrupting functions and them to the observation composite type
    rngcorrupt = Xoshiro(corruptseed)
    rngtrop = Xoshiro(tropseed)

    # construct CjlObservation object    
    observation = CjlObservation{Float64}(msname,data3dresandperm,flag3dresandperm,flagrow,antenna1,antenna2,uvw,times,exposure,scanno,numchan,chanfreqvec,
    chanwidth,phasedir,pos,stationinfo,tropwetonly,correff,tropattenuate,tropskynoise,tropmeandelays,tropturbulence,polframe,polmode,
    ptginterval,ptgscale,ptgmode,stationgainsmode,bandpassfile,rngcorrupt,rngtrop)
=#
