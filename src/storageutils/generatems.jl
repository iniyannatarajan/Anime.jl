export msfromconfig, msfromuvfits

using StatsBase: mode

"""
    makecasaanttable(stations::String, casaanttemplate::String; delim::String=",", ignorerepeated::Bool=false)

Generate a CASA antenna table from CSV station info file in the current working directory.
"""
function makecasaanttable(stations::String, casaanttemplate::String; delim::String=",", ignorerepeated::Bool=false)
    # read in the stations CSV file
    df = CSV.read(stations, DataFrame; delim=delim, ignorerepeated=ignorerepeated)
    
    stationtable = "ANTENNA_$(split(basename(stations), '.')[1])"

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
    addweightcols(msname::String, mode::String, sigmaspec::Bool, weightspec::Bool)

Add WEIGHT_SPECTRUM and SIGMA_SPECTRUM columns to the MS
"""
function addweightcols(msname::String, mode::String, sigmaspec::Bool, weightspec::Bool)
    # TODO get mode as arg and if it is "uvfits", delete wtspec and regenerate
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
    msfromuvfits(uvfits::String, msname::String, mscreationmode::String, stations::String; delim::String=",", ignorerepeated::Bool=false)

Generate MS from existing UVFITS file
"""
function msfromuvfits(uvfits::String, msname::String, mscreationmode::String, stations::String; delim::String=",", ignorerepeated::Bool=false)
    # convert uvfits to ms
    importuvfits(fitsfile=uvfits, vis=msname)

    # compare ANTENNA table in ms with stations file
    df = CSV.read(stations, DataFrame; delim=delim, ignorerepeated=ignorerepeated)

    tb.open("$(msname)::ANTENNA")
    msstations = string.(tb.getcol("STATION"))
    tb.close()

    if !(issetequal(msstations, df.station))
        error("Station info file $(stations) and MS ANTENNA table do not match 🤷")
    end
   
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
    addweightcols(msname, mscreationmode, true, true)
end

"""
    msfromconfig(msname::String, mscreationmode::String, stations::String, casaanttemplate::String, spw_centrefreq::Array{Float64, 1}, 
    spw_bw::Array{Float64, 1}, spw_channels::Array{Int64, 1}, sourcedict::Dict{String, Any}, starttime::String, exposure::Float64, scans::Int64,
    scanlengths::Array{Float64, 1}, scanlag::Float64; autocorr::Bool=false, telescopename::String="VLBA", feed::String="perfect R L",
    shadowlimit::Float64=1e-6, elevationlimit::String="10deg", stokes::String="RR RL LR LL", delim::String=",", ignorerepeated::Bool=false)

Generate measurement set from scratch from input observation parameters
"""
function msfromconfig(msname::String, mscreationmode::String, stations::String, casaanttemplate::String, spw_centrefreq::Array{Float64, 1}, 
    spw_bw::Array{Float64, 1}, spw_channels::Array{Int64, 1}, sourcedict::Dict{String, Any}, starttime::String, exposure::Float64, scans::Int64,
    scanlengths::Array{Float64, 1}, scanlag::Float64; autocorr::Bool=false, telescopename::String="VLBA", feed::String="perfect R L",
    shadowlimit::Float64=1e-6, elevationlimit::String="10deg", stokes::String="RR RL LR LL", delim::String=",", ignorerepeated::Bool=false)
    # open new MS
    sm.open(msname)

    # set autocorr
    Bool(autocorr) ? sm.setauto(autocorrwt=1.0) : sm.setauto(autocorrwt=0.0)
    
    # set antenna config -- create a new CASA ANTENNA table if the station info is in a CSV file
    #stationtable = ""
    if isfile(stations) && isdir(casaanttemplate)
        @info("Creating new ANTENNA table from CSV station info file...")
	    stationtable = makecasaanttable(stations, casaanttemplate; delim=delim, ignorerepeated=ignorerepeated)
    else
        error("Either $(stations) or $(casaanttemplate) does not exist 🤷")
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

    addweightcols(msname, mscreationmode, true, true) # add weight columns

end