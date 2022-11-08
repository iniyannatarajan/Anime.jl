using CSV
using DataFrames

# import python libraries
simulator = pyimport("casatools" => "simulator")
table = pyimport("casatools" => "table")
measures = pyimport("casatools" => "measures")

# instantiate casa python classes
sm = simulator()
tb = table()
me = measures()

function mkCasaAntTable(stationasciifile::String, delim::String, ignorerepeated::Bool, casaanttemplate::String)
    """
    Generate a CASA antenna table from CSV station info file in the current working directory.
    """
    # read in the ASCII file
    df = CSV.read(stationasciifile, DataFrame; delim=delim, ignorerepeated=ignorerepeated)
    
    stationtable = "ANTENNA_$(split(basename(stationasciifile), '.')[1])"

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

function genms(yamlconf::Dict, casaanttemplate::String)
    """
    Main function to generate MS.
    """
    # open new MS
    sm.open(yamlconf["msname"])

    # set autocorr
    Bool(yamlconf["manual"]["autocorr"]) ? sm.setauto(autocorrwt=1.0) : sm.setauto(autocorrwt=0.0)
    
    # set antenna config -- optionally create a new CASA ANTENNA table if the station info is in a CSV file
    stationtable = ""
    if isfile(yamlconf["manual"]["stations"])
	# check if template is specified
	casaanttemplate == nothing && error("$(yamlconf["manual"]["stations"]) is a CSV file but template CASA ANTENNA table not specified!")
	@info("Creating new ANTENNA table from CSV station info file...")
	stationtable = mkCasaAntTable(yamlconf["manual"]["stations"], " ", true, casaanttemplate)
	@info("ANTENNA table $(stationtable) created")
    elseif isdir(yamlconf["manual"]["stations"])
	stationtable = yamlconf["manual"]["stations"]
    else
        error("$(yamlconf["manual"]["stations"]) does not exist!")
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
    obspos = me.observatory(yamlconf["telescopename"])
    me.doframe(obspos)

    # NB: For any Python array with strings, convert elements to Julia strings and cast the entire Vector to PyList()
    sm.setconfig(telescopename=yamlconf["telescopename"], x=x, y=y, z=z,
		 dishdiameter=dish_diameter, mount=PyList(string.(station_mount)), antname=PyList(string.(station_ind)),
		 padname=PyList(string.(station_code)), coordsystem=coords, referencelocation=obspos)

    # set feed
    sm.setfeed(mode=yamlconf["manual"]["feed"])

    # set limits for data flagging
    sm.setlimits(shadowlimit=yamlconf["manual"]["shadowlimit"], elevationlimit=yamlconf["manual"]["elevationlimit"])

    # set spectral windows
    for (key, val) in yamlconf["manual"]["spwname"]
	startfreq = val["centrefreq"] - val["bandwidth"]-2.0
	chanwidth = val["bandwidth"]/val["channels"]
	sm.setspwindow(spwname=key, freq="$(startfreq)GHz",
		deltafreq="$(chanwidth)GHz", freqresolution="$(chanwidth)GHz",
		nchannels=val["channels"], stokes=yamlconf["manual"]["stokes"])
    end

    # set obs fields
    for (key, val) in yamlconf["manual"]["source"]
	sm.setfield(sourcename=key, sourcedirection=me.direction(rf=val["epoch"], v0=val["RA"], v1=val["Dec"]))
    end

    # set obs times
    obs_starttime = split(yamlconf["manual"]["starttime"], ",")
    referencetime = me.epoch(obs_starttime...)
    me.doframe(referencetime)
    sm.settimes(integrationtime=yamlconf["manual"]["inttime"], usehourangle=false, referencetime=referencetime)

    # observe
    Int64(yamlconf["manual"]["scans"]) != length(yamlconf["manual"]["scanlengths"]) && error("Number of scans and length(scanlengths_s) do not match!")
    Int64(yamlconf["manual"]["scans"]) != length(yamlconf["manual"]["scanlags"])+1 && error("Number of scans and length(scanlaglist)+1 do not match!")

    stoptime = 0
    for ii in 1:Int64(yamlconf["manual"]["scans"])
	starttime = ii==1 ? 0 : stoptime+yamlconf["manual"]["scanlags"][ii-1]
	stoptime = starttime + yamlconf["manual"]["scanlengths"][ii]

	for (key, val) in yamlconf["manual"]["spwname"]
	    sm.observe(sourcename=collect(keys(yamlconf["manual"]["source"]))[1], spwname=key, starttime="$(starttime)s", stoptime="$(stoptime)s")
	end
    end

    # clean up
    sm.close()
    @info("$(yamlconf["msname"]) successfully created.")
end
