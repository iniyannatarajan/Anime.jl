export stationgains

function generatetimeseries(mode::String, location::ComplexF32, scale::Float64, nsamples::Int64)
    """
    Generate complex-valued wiener series
    """
    series = zeros(ComplexF32, nsamples)
    if mode == "wiener"	    
        sqrtnsamples = sqrt(nsamples)
        series[1] = location + (scale*randn(ComplexF32)/sqrtnsamples)
        for ii in 2:nsamples
            series[ii] = series[ii-1] + (scale*randn(ComplexF32)/sqrtnsamples)
        end
    elseif mode == "gaussian"
	series = location + scale*randn(ComplexF32, nsamples)
    end
    return series
end

function stationgains(obs::CjlObservation)
    """
    Add time-vriable station gains to visibilities
    """
    # get element type to be used
    elemtype = typeof(obs.data[1][1])

    # open h5 file for writing
    fid = h5open(obs.yamlconf["hdf5corruptions"], "r+")
    g = create_group(fid, "stationgains")
    attributes(g)["desc"] = "Numerical values of time-variable per station G-Jones terms applied to data"
    attributes(g)["format"] = "for each scan, a 4d array of 2 x 2 x ntime x nant is stored"

    # get unique scan numbers
    uniqscans = unique(obs.scanno)

    # loop through each station to create a vector of 2x2 G-Jones terms evolving over time
    gjonesdict = Dict() # create empty dict
    for scan in uniqscans
        # compute ideal ntimes per scan
        actualtscanvec = unique(getindex(obs.times, findall(obs.scanno.==scan)))
	idealtscanvec = collect(first(actualtscanvec):obs.exposure:last(actualtscanvec))
	idealtscanveclen = length(idealtscanvec)

	# create 4d array to hold G-Jones terms per time per station
	gjonesmatrices = zeros(elemtype, 2, 2, idealtscanveclen, size(obs.stationinfo)[1]) # 2 x 2 x ntimes x nant

	for ant in 1:size(obs.stationinfo)[1]
	   gjonesmatrices[1, 1, :, ant] = generatetimeseries(obs.yamlconf["stationgains"]["mode"], obs.stationinfo.g_pol1_loc[ant], obs.stationinfo.g_pol1_scale[ant], idealtscanveclen)
	   gjonesmatrices[2, 2, :, ant] = generatetimeseries(obs.yamlconf["stationgains"]["mode"], obs.stationinfo.g_pol2_loc[ant], obs.stationinfo.g_pol2_scale[ant], idealtscanveclen)
	end

	# add to gjonesdict
	gjonesdict[scan] = gjonesmatrices

        # loop over time/row and apply gjones terms corresponding to each baseline
	for idealtimeindex in 1:idealtscanveclen
	    if idealtscanvec[idealtimeindex] in actualtscanvec
		#for chan in 1:obs.numchan
		    #obs.data[:,:,chan,idealtimeindex] = gjonesdict[scan][:,:,idealtimeindex,obs.antenna1[idealtimeindex]+1]*obs.data[:,:,chan,idealtimeindex]*adjoint(gjonesdict[scan][:,:,idealtimeindex,obs.antenna2[idealtimeindex]+1])
		    obs.data[:,:,:,idealtimeindex] = gjonesdict[scan][:,:,idealtimeindex,obs.antenna1[idealtimeindex]+1].*obs.data[:,:,:,idealtimeindex].*adjoint(gjonesdict[scan][:,:,idealtimeindex,obs.antenna2[idealtimeindex]+1])
	        #end
	    else
		continue
	    end
        end

	# write to h5 file
	g["scan $(scan)"] = gjonesmatrices
    end

    # add datatype attribute
    attributes(g)["datatype"] = string(typeof(read(g[keys(g)[1]])))

    # close h5 file
    close(fid)

    @info("Compute and apply station gains... 🙆")
end
