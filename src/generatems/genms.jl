using CSV
using DataFrames
using JSON3

# import python libraries
simulator = pyimport("casatools" => "simulator")
table = pyimport("casatools" => "table")
measures = pyimport("casatools" => "measures")

# instantiate casa python classes
sm = simulator()
tb = table()
me = measures()

function mkCasaAntTable(stationasciifile::String, delim::String, ignorerepeated::Bool, casatemplate::String)
    """
    Generate a CASA antenna table from CSV station info file in the current working directory.
    """
    # read in the ASCII file
    df = CSV.read(stationasciifile, DataFrame; delim=delim, ignorerepeated=ignorerepeated)
    
    stationtable = "ANTENNA_$(split(basename(stationasciifile), '.')[1])"

    # read template table and copy
    tb.open(casatemplate)
    tb.copy(stationtable, deep=true, valuecopy=true)
    tb.close()
    tb.clearlocks()

    # populate new table with values from the DataFrame
    tb.open(stationtable, nomodify=false)
    for ii in 1:length(df.station)-1
        tb.copyrows(stationtable, nrow=1)
    end
    tb.putcol("STATION", PyList(string.(df.station)))
    tb.putcol("NAME", PyList(string.(1:length(df.station)))) # assign numbers as indices/names for the stations
    tb.putcol("MOUNT", PyList(string.(df.mount)))
    tb.putcol("DISH_DIAMETER", PyList(df.dish_diameter))
    for jj in 1:length(df.station)
	tb.putcell("POSITION", jj-1, PyList(parse.(Float64, split.(df.xzy_position_m,',')[jj])))
    end
    tb.close()
    tb.clearlocks()

    return stationtable
end

function genms(jsonpars::JSON3.Object, templatetable::String, clobber::Bool)
    """
    Main function to generate MS.
    """
    # check if MS exists and optionally delete it
    isdir(jsonpars.msname) ? (clobber || error("$(jsonpars.msname) exists! Not overwriting.")) : run(`rm -rf $(jsonpars.msname)`)

    # open new MS
    sm.open(jsonpars.msname)

    # set autocorr
    Bool(jsonpars.toggle_autocorr) ? sm.setauto(autocorrwt=1.0) : sm.setauto(autocorrwt=0.0)
    
    # set antenna config -- optionally create a new CASA ANTENNA table if the station info is in a CSV file
    stationtable = ""
    if isfile(jsonpars.stationtable)
	# check if template is specified
	templatetable == nothing && error("$(jsonpars.stationtable) is a CSV file but template CASA ANTENNA table not specified!")
	@info("Creating new ANTENNA table from CSV station info file...")
	stationtable = mkCasaAntTable(jsonpars.stationtable, " ", true, templatetable)
	@info("ANTENNA table $(stationtable) created")
    elseif isdir(jsonpars.stationtable)
	stationtable = jsonpars.stationtable
    else
        error("$(jsonpars.stationtable) does not exist!")
    end

    # station locations are contained in a casa antenna table
    tb.open(stationtable)
    station_code = tb.getcol("STATION")
    station_ind = tb.getcol("NAME")
    dish_diameter = tb.getcol("DISH_DIAMETER")
    station_mount = tb.getcol("MOUNT")
    x, y, z = tb.getcol("POSITION")
    tb.close()
    tb.clearlocks()

    coords = "global" # CASA antenna tables use ITRF by default
    obspos = me.observatory(jsonpars.telname)
    me.doframe(obspos)

    # NB: For any Python array with strings, convert elements to Julia strings and cast the entire Vector to PyList()
    sm.setconfig(telescopename=jsonpars.telname, x=x, y=y, z=z,
		 dishdiameter=dish_diameter, mount=PyList(string.(station_mount)), antname=PyList(string.(station_ind)),
		 padname=PyList(string.(station_code)), coordsystem=coords, referencelocation=obspos)

    # set feed
    sm.setfeed(mode=jsonpars.feed)

    # set limits for data flagging
    sm.setlimits(shadowlimit=jsonpars.shadowlimit_frac, elevationlimit=jsonpars.elevationlimit)

    # set spectral windows
    length(jsonpars.centrefreq_ghz) != length(jsonpars.spwnames) && error("Check if centre frequencies are defined for all spectral windows or vice versa!")
    length(jsonpars.bandwidth_ghz) != length(jsonpars.spwnames) && error("Check if bandwidths are defined for all spectral windows or vice versa!")

    for ii in 1:length(jsonpars.spwnames)
	startfreq = jsonpars.centrefreq_ghz[ii] - jsonpars.bandwidth_ghz[ii]/2
	deltafreq = jsonpars.bandwidth_ghz[ii]/jsonpars.nchan
	sm.setspwindow(spwname=jsonpars.spwnames[ii], freq="$(startfreq)GHz",
		deltafreq="$(deltafreq)GHz", freqresolution="$(deltafreq)GHz",
                nchannels=jsonpars.nchan, stokes=jsonpars.poltype)
    end

    # set obs fields
    for (src, val) in jsonpars.sources
        sm.setfield(sourcename=src, sourcedirection=me.direction(rf=val.epoch, v0=val.ra, v1=val.dec))
    end

    # set obs times
    obs_starttime = split(jsonpars.starttime, ",")
    referencetime = me.epoch(obs_starttime...)
    me.doframe(referencetime)
    sm.settimes(integrationtime=jsonpars.inttime_s, usehourangle=false, referencetime=referencetime)

    # observe
    length(jsonpars.scans) != length(jsonpars.scanlengths_s) && error("scanlist and scanlengthlist do not match!")
    length(jsonpars.scans) != length(jsonpars.scanlags_s)+1 && error("scanlist and scanlaglist+1 do not match!")

    stoptime = 0
    for ii in 1:length(jsonpars.scans)
	starttime = ii==1 ? 0 : stoptime+jsonpars.scanlags_s[ii-1]
	stoptime = starttime + jsonpars.scanlengths_s[ii]

	for jj in 1:length(jsonpars.spwnames)
    	    sm.observe(sourcename=jsonpars.scans[ii], spwname=jsonpars.spwnames[jj], starttime="$(starttime)s", stoptime="$(stoptime)s")
	end
    end

    # clean up
    sm.close()
    @info("$(jsonpars.msname) successfully created.")
end
