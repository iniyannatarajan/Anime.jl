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
	series = location .+ scale*randn(ComplexF32, nsamples)
    end
    return series
end

function stationgains(obs::CjlObservation)
    """
    Add time-variable station gains to visibilities
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
    row = 1 # variable to index obs.data array
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
	    currenttime = idealtscanvec[idealtimeindex]
	    if currenttime in actualtscanvec
		# read all baselines present in a given time
		ant1vec = getindex(obs.antenna1, findall(obs.times.==currenttime))
		ant2vec = getindex(obs.antenna2, findall(obs.times.==currenttime))
		for (ant1,ant2) in zip(ant1vec, ant2vec)
		    for chan in 1:obs.numchan
			#@info(row, scan, idealtimeindex, chan, ant1+1, ant2+1)
                        #@info(gjonesdict[scan][:,:,idealtimeindex,ant1+1], obs.data[:,:,chan,row], adjoint(gjonesdict[scan][:,:,idealtimeindex,ant2+1]))
		        obs.data[:,:,chan,row] = gjonesdict[scan][:,:,idealtimeindex,ant1+1]*obs.data[:,:,chan,row]*adjoint(gjonesdict[scan][:,:,idealtimeindex,ant2+1])
                        #@info(obs.data[:,:,chan,row])
		    end
		    row += 1 # increment obs.data last index i.e. row number
	        end
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
