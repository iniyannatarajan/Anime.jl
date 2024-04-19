export makecasaantennatable, addweightcolumns, createmsfromuvfits, createuvfitsfromms, createmsfromconfig, writems

using StatsBase: mode

"""
    makecasaantennatable(df::DataFrame, casaanttemplate::String; delim::String=",", ignorerepeated::Bool=false)

Create CASA antenna table from the DataFrame df containing station information using `casaanttemplate` as template. Requires `casatools`.
"""
function makecasaantennatable(df::DataFrame, casaanttemplate::String; delim::String=",", ignorerepeated::Bool=false)    
    #stationtable = "ANTENNA_$(split(basename(stations), '.')[1])"
    stationtable = "ANTENNA_temp"

    # read template table and copy
    tb.open(casaanttemplate)
    tb.copy(stationtable, deep=true, valuecopy=true)
    tb.close()
    tb.clearlocks()

    # populate new table with values from the DataFrame
    tb.open(stationtable, nomodify=false)
    for ii in 1:length(df.station)-1
        tb.copyrows(stationtable, nrow=1)
    end
    for ii in 1:length(df.station)
        tb.putcol("STATION", String(df.station[ii]), startrow=ii-1, nrow=1)
        tb.putcol("NAME", String(df.station[ii]), startrow=ii-1, nrow=1)
        tb.putcol("MOUNT", String(df.mount[ii]), startrow=ii-1, nrow=1)
    end
    tb.putcol("DISH_DIAMETER", PyList(df.dishdiameter_m))
    #for jj in 1:length(df.station)
    #	tb.putcell("POSITION", jj-1, PyList(parse.(Float64, split.(df.xyzpos_m,',')[jj])))
    #end
    index = 0
    for (x,y,z) in zip(df.x_m,df.y_m,df.z_m)
        tb.putcell("POSITION", index, PyList([x,y,z]))
        index += 1
    end
    tb.close()
    tb.clearlocks()

    return stationtable
end

"""
    addweightcolumns(msname::String, mode::String, sigmaspec::Bool, weightspec::Bool)

Add ```WEIGHT_SPECTRUM``` and ```SIGMA_SPECTRUM``` columns to existing MS. Requires `casatools` to be installed.
"""
function addweightcolumns(msname::String, mode::String, sigmaspec::Bool, weightspec::Bool)
    # get quantities to define array shapes
    tb.open("$(msname)::SPECTRAL_WINDOW")
    numchan = pyconvert(Int64, tb.getcol("NUM_CHAN")[0])
    tb.close()
    
    tb.open("$(msname)::POLARIZATION")
    numcorr = pyconvert(Int64, tb.getcol("NUM_CORR")[0])
    tb.close()

    shape = PyList([numchan, numcorr])

    tb.open(msname, nomodify=false)
    if sigmaspec
        coldesc = pydict(Dict("valueType" => "float", "dataManagerType" => "TiledShapeStMan", "dataManagerGroup" => "TiledSigmaSpectrum", "option" => 0, "maxlen" => 0, "comment" => "Estimated rms noise for each data point", "ndim" => 2, "_c_order" => true, "keywords" => pydict(Dict()), "shape" => shape))
        tabdesc = pydict(Dict("SIGMA_SPECTRUM" => coldesc))
	tb.addcols(tabdesc)
	sigsp = tb.getcol("SIGMA_SPECTRUM")
	PyArray(sigsp)[:] .= 1.0
	tb.putcol("SIGMA_SPECTRUM", sigsp)
    end

    if mode == "uvfits"
        tb.removecols("WEIGHT_SPECTRUM") # remove automatically generated weight spectrum column to avoid dimension mismatch
    end

    if weightspec
	coldesc = pydict(Dict("valueType" => "float", "dataManagerType" => "TiledShapeStMan", "dataManagerGroup" => "TiledWgtSpectrum", "option" => 0, "maxlen" => 0, "comment" => 
			      "Weight for each data point", "ndim" => 2, "_c_order" => true, "keywords" => pydict(Dict()), "shape" => shape))
	tabdesc = pydict(Dict("WEIGHT_SPECTRUM" => coldesc))
	tb.addcols(tabdesc)
	wtsp = tb.getcol("WEIGHT_SPECTRUM")
	PyArray(wtsp)[:] .= 1.0
	tb.putcol("WEIGHT_SPECTRUM", wtsp)
    end
    tb.close()
end

#="""
    setup_polarization(msname::String, stationtable::String)

Set up all polarization subtables in the MS based on the station feed types
"""
function setup_polarization(msname::String, stationtable::String)
    # get feed types from the stationtable
    #df = CSV.read(stationtable, DataFrame; delim=delim, ignorerepeated=ignorerepeated)
    #feedtypes = df["feed"]

    # modify DATA_DESCRIPTION, POLARIZATION, FEED subtables

end=#

#="""
    msfromvex()

Extract config parameters from VEX schedule and call msfromconfig() to generate MS
"""
function msfromvex()
end=#

"""
    createmsfromuvfits(uvfits::String, msname::String, mscreationmode::String; overwrite::Bool=true)

Convert `uvfits` to MS. Requires `casatools` and `casatasks`.
"""
function createmsfromuvfits(uvfits::String, msname::String, mscreationmode::String; overwrite::Bool=true)

    # overwrite existing MS
    if isdir(msname)
        if !overwrite
            @error("$msname exists and overwrite set to $overwrite; not converting uvfits to ms 🤷")
            exit()
        else
            @info("Removing existing $msname...")
            rm(msname; recursive=true)
       end
    end

    # convert uvfits to ms
    importuvfits(fitsfile=uvfits, vis=msname)

    #= # compare ANTENNA table in ms with stations file
    df = CSV.read(stations, DataFrame; delim=delim, ignorerepeated=ignorerepeated)

    tb.open("$(msname)::ANTENNA")
    msstations = string.(tb.getcol("STATION"))
    tb.close()

    if !(issetequal(msstations, df.station))
        error("Station info file $(stations) and MS ANTENNA table do not match 🤷")
    end=#
   
    # setup polarization info
    #setup_polarization(msname, stations)

    # replace EXPOSURE column with the mode of EXPOSURE column in the ms to avoid slight differences
    # in recorded exposure times in uvfits output by eht-imaging (the most likely origin of a uvfits file)
    tb.open(msname, nomodify=false)
    exparr = pyconvert(PyArray, tb.getcol("EXPOSURE"))
    @info("Replacing potentially inconsistent values in EXPOSURE and INTERVAL columns with mode(EXPOSURE): $(mode(exparr)) s")
    exposurevec = mode(exparr) .+ zeros(Float64, pyconvert(Int64, tb.nrows())) # use mode instead of mean or median to get the "most correct" exposure value
    tb.putcol("EXPOSURE", PyList(exposurevec))
    tb.putcol("INTERVAL", PyList(exposurevec))
    tb.close()

    # WEIGHT_SPECTRUM is added by importuvfits; add only SIGMA_SPECTRUM manually
    addweightcolumns(msname, mscreationmode, true, true)
end

"""
    createuvfitsfromms(msname::String, uvfits::String, datacolumn::String; field::String="", spw::String="", antenna::String="",
    timerange::String="", overwrite::Bool=true)

Create UVFITS from existing MS. Requires `casatasks`.
"""
function createuvfitsfromms(msname::String, uvfits::String, datacolumn::String; field::String="", spw::String="", antenna::String="",
    timerange::String="", overwrite::Bool=true)
    @info("Creating $uvfits from $msname...")

    if !overwrite && isfile(uvfits)
        @error("$uvfits exists but overwrite=$overwrite; not exporting to uvfits 🤷")
        exit()
    end
    # convert ms to uvfits
    exportuvfits(vis=msname, fitsfile=uvfits, datacolumn=datacolumn, field=field, spw=spw, antenna=antenna, timerange=timerange, overwrite=overwrite)
end

"""
    createmsfromconfig(msname::String, mscreationmode::String, casaanttemplate::String, stationinfo::DataFrame, spw_centrefreq::Array{Float64, 1}, 
    spw_bw::Array{Float64, 1}, spw_channels::Array{Int64, 1}, sourcedict::Dict{String, Any}, starttime::String, exposure::Float64, scans::Int64,
    scanlengths::Array{Float64, 1}, scanlag::Float64; autocorr::Bool=false, telescopename::String="VLBA", feed::String="perfect R L",
    shadowlimit::Float64=1e-6, elevationlimit::String="10deg", stokes::String="RR RL LR LL", delim::String=",", ignorerepeated::Bool=false)

Create Measurement Set from observation parameters and DataFrame containing station information. Requires `casatools`.
"""
function createmsfromconfig(msname::String, mscreationmode::String, casaanttemplate::String, stationinfo::DataFrame, spw_centrefreq::Array{Float64, 1}, 
    spw_bw::Array{Float64, 1}, spw_channels::Array{Int64, 1}, sourcedict::Dict{String, Any}, starttime::String, exposure::Float64, scans::Int64,
    scanlengths::Array{Float64, 1}, scanlag::Float64; autocorr::Bool=false, telescopename::String="VLBA", feed::String="perfect R L",
    shadowlimit::Float64=1e-6, elevationlimit::String="10deg", stokes::String="RR RL LR LL", delim::String=",", ignorerepeated::Bool=false)
    # open new MS
    sm.open(msname)

    # set autocorr
    Bool(autocorr) ? sm.setauto(autocorrwt=1.0) : sm.setauto(autocorrwt=0.0)
    
    # set antenna config -- create a new CASA ANTENNA table if the station info is in a CSV file
    if isdir(casaanttemplate)
        @info("Creating new ANTENNA table using stationinfo DataFrame...")
	    stationtable = makecasaantennatable(stationinfo, casaanttemplate; delim=delim, ignorerepeated=ignorerepeated)
    else
        error("The CASA ANTENNA template table $(casaanttemplate) does not exist 🤷")
    end

    # station locations are contained in a casa antenna table
    tb.open(stationtable)
    station_code = tb.getcol("STATION")
    dish_diameter = tb.getcol("DISH_DIAMETER")
    station_mount = tb.getcol("MOUNT")
    x, y, z = tb.getcol("POSITION")
    tb.close()
    tb.clearlocks()

    coords = "global" # CASA antenna tables use ITRF by default
    obspos = me.observatory(telescopename)
    me.doframe(obspos)

    # NB: For any Python array with strings, convert elements to Julia strings and cast the entire Vector to PyList()
    sm.setconfig(telescopename=telescopename, x=x, y=y, z=z,
		 dishdiameter=dish_diameter, mount=PyList(string.(station_mount)), antname=PyList(string.(station_code)),
		 padname=PyList(string.(station_code)), coordsystem=coords, referencelocation=obspos)

    # set feed
    sm.setfeed(mode=feed)

    # set limits for data flagging
    sm.setlimits(shadowlimit=shadowlimit, elevationlimit=elevationlimit)

    # set spectral windows -- NB: we are not handling multiple spws as of now
    #=for (key, val) in yamlconf["manual"]["spwname"]
	startfreq = val["centrefreq"] - val["bandwidth"]-2.0
	chanwidth = val["bandwidth"]/val["channels"]
	sm.setspwindow(spwname=key, freq="$(startfreq)GHz",
		deltafreq="$(chanwidth)GHz", freqresolution="$(chanwidth)GHz",
		nchannels=val["channels"], stokes=yamlconf["manual"]["stokes"])
    end=#
    for ind in 1:length(spw_centrefreq)
        chanwidth = spw_bw[ind]/spw_channels[ind]
        startfreq = spw_centrefreq[ind] - spw_bw[ind]/2. + chanwidth/2.
	sm.setspwindow(spwname="$(Int(spw_centrefreq[ind]/1e9))GHz", freq="$(startfreq)Hz",
        	deltafreq="$(chanwidth)Hz", freqresolution="$(chanwidth)Hz",
		nchannels=spw_channels[ind], stokes=stokes)
    end

    # set obs fields
    for (key, val) in sourcedict
	sm.setfield(sourcename=key, sourcedirection=me.direction(rf=val["epoch"], v0=val["RA"], v1=val["Dec"]))
    end

    # set obs times
    obs_starttime = split(starttime, ",")
    referencetime = me.epoch(obs_starttime...)
    me.doframe(referencetime)
    sm.settimes(integrationtime=exposure, usehourangle=false, referencetime=referencetime)

    # observe
    scans != length(scanlengths) && error("No. of scans and no. of scanlengths do not match 🤷")

    stoptime = 0
    for ii in 1:scans
	    starttime = ii==1 ? 0 : stoptime+scanlag
	    stoptime = starttime + scanlengths[ii] #+ yamlconf["exposure"] # one more exposure added here at the end of the scan

	    # NB: We are not handling multiple spectral windows as of now
	    #=for (key, val) in yamlconf["manual"]["spwname"]
	    sm.observe(sourcename=collect(keys(yamlconf["manual"]["source"]))[1], spwname=key, starttime="$(starttime)s", stoptime="$(stoptime)s")
	    end=#
	    for ind in 1:length(spw_centrefreq)
            sm.observe(sourcename=collect(keys(sourcedict))[1], spwname="$(Int(spw_centrefreq[ind]/1e9))GHz", starttime="$(starttime)s", stoptime="$(stoptime)s")
	    end
    end

    # clean up
    sm.close()

    addweightcolumns(msname, mscreationmode, true, true) # add weight columns

end

"""
    createmsfromconfig(obsconfig, stationinfo::DataFrame)

Alias with fewer arguments for use in pipelines.
"""
function createmsfromconfig(obsconfig::Dict, stationinfo::DataFrame)
    createmsfromconfig(obsconfig["msname"], obsconfig["mode"], obsconfig["casaanttemplate"], stationinfo, obsconfig["spw"]["centrefreq"], obsconfig["spw"]["bandwidth"],
    obsconfig["spw"]["channels"], obsconfig["source"], obsconfig["starttime"], obsconfig["exposure"], obsconfig["scans"], obsconfig["scanlengths"], obsconfig["scanlag"];
    autocorr=obsconfig["autocorr"], telescopename=obsconfig["telescopename"], feed=obsconfig["feed"], shadowlimit=obsconfig["shadowlimit"],
    elevationlimit=obsconfig["elevationlimit"], stokes=obsconfig["stokes"], delim=",", ignorerepeated=false)
end

"""
    computeweights!(totalrmsspec::Array{Float32, 4}, totalwtspec::Array{Float32, 4}; h5file::String="")

Compute total rms (sigma) values and inverse-squared visibility weights from thermal+sky noise terms stored in `h5file`.
"""
function computeweights!(totalrmsspec::Array{Float32, 4}, totalwtspec::Array{Float32, 4}; h5file::String="")

    fid = h5open(h5file, "r")
    if haskey(fid, "thermalnoise") && haskey(fid["thermalnoise"], "thermalnoiserms")
        totalrmsspec[:,:,:,:] = read(fid["thermalnoise"]["thermalnoiserms"]).^2
    end
    if haskey(fid, "troposphere") && haskey(fid["troposphere"], "skynoiserms")
        totalrmsspec[:,:,:,:] = totalrmsspec[:,:,:,:] + read(fid["troposphere"]["skynoiserms"]).^2
    end

    # compute sigma_spectrum
    totalrmsspec[:,:,:,:] = sqrt.(totalrmsspec[:,:,:,:])

    # compute weight_spectrum
    totalwtspec[:,:,:,:] = 1 ./(totalrmsspec[:,:,:,:].^2)

    return totalrmsspec, totalwtspec
end

"""
    writems(ms::MeasurementSet; h5file::String="")

Write data to MS. Optionally, add weight columns from `h5file`.
"""
function writems(ms::MeasurementSet; h5file::String="")
    # replace NaNs with zeros
    ms.data[isnan.(ms.data)] .= 0.0+0.0*im

    # compute weight columns
    #@warn("WEIGHT_SPECTRUM and SIGMA_SPECTRUM are not filled in properly at the moment, while the code is being optimised.")
    totalrmsspec = zeros(Float32, 2, 2, ms.numchan, size(ms.data)[4]) # noise rms zero by default
    totalwtspec = ones(Float32, 2, 2, ms.numchan, size(ms.data)[4]) # weights unity by default

    if isfile(h5file)
        @info("Populating weight and sigma spectrum arrays...")
        computeweights!(totalrmsspec, totalwtspec, h5file=h5file)
    end

    # convert the sigma_spec and weight_spec arrays to the format required by Casacore.jl
    revtotalrmsspec3dres = permutedims(totalrmsspec, (2,1,3,4))
    revtotalrmsspec3d = reshape(revtotalrmsspec3dres, 4, size(revtotalrmsspec3dres)[3], :)
    revtotalrmsspec = [Matrix{Float32}(revtotalrmsspec3d[:,:,i]) for i in 1:size(revtotalrmsspec3d)[3]]

    revtotalwtspec3dres = permutedims(totalwtspec, (2,1,3,4))
    revtotalwtspec3d = reshape(revtotalwtspec3dres, 4, size(revtotalwtspec3dres)[3], :)
    revtotalwtspec = [Matrix{Float32}(revtotalwtspec3d[:,:,i]) for i in 1:size(revtotalwtspec3d)[3]]

    # compute channel averaged weight and sigma column values
    revtotalrms = Vector{Vector{Float32}}()
    revtotalwt = Vector{Vector{Float32}}()
    for row in 1:size(ms.data)[4] # nrows
        tmpvec = dropdims(mean(revtotalrmsspec[row], dims=2), dims=2)
        push!(revtotalrms, tmpvec) # average the sigma values
        push!(revtotalwt, 1 ./tmpvec.^2)
    end

    # reshape data array for Casacore.jl
    revdata3dres = permutedims(ms.data, (2,1,3,4))
    revdata3d = reshape(revdata3dres, 4, size(revdata3dres)[3], :)
    revdata = [Matrix{ComplexF32}(revdata3d[:,:,i]) for i in 1:size(revdata3d)[3]]

    # when all the corruptions have been applied, write the above columns back to disk
    table = CCTable(ms.name, CCTables.Update)
    table[:DATA] = revdata
    table[:SIGMA_SPECTRUM] = revtotalrmsspec
	table[:WEIGHT_SPECTRUM] = revtotalwtspec
    table[:SIGMA] = revtotalrms
	table[:WEIGHT] = revtotalwt

    @info("Write arrays to disk... 🙆")
end