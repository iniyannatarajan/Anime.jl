export instrumentalpol

include(joinpath("util.jl"))

function instrumentalpol(obs::CjlObservation)
    """
    compute and apply polarization leakage
    """
    # get element type to be used for the Jones matrices
    elemtype = typeof(obs.data[1][1])

    # get unique times
    uniqtimes = unique(obs.times)

    # open h5 file for writing
    fid = h5open(obs.yamlconf["hdf5corruptions"], "r+")
    g = create_group(fid, "polarization")
    attributes(g)["desc"] = "Numerical values of instrumental polarization matrices applied to data (Dterms and fullpolproduct)"
    attributes(g)["dims"] = "2 x 2 x nchannels x nant" #for each scan, a 4d array of 2 x 2 x nchan x nant is stored

    # compute necessary quantities
    elevationmatrix = elevationangle(obs)
    parallacticanglematrix = parallacticangle(obs)

    # compute D-terms
    if obs.yamlconf["instrumentalpol"]["visibilityframe"] != "antenna"
	# compute in sky frame
	if obs.yamlconf["instrumentalpol"]["visibilityframe"] != "sky"
	    @warn("Visibility frame not recognised! Computing visibilities in sky frame...")
	end
        
	# D-terms -- perform twice the feed angle rotation
	djonesmatrices = ones(elemtype, 2, 2, obs.numchan, size(obs.stationinfo)[1]) # 2 x 2 x nchan x nant
	polrotmatrices = ones(elemtype, 2, 2, obs.numchan, size(uniqtimes)[1], size(obs.stationinfo)[1])

        for ant in 1:size(obs.stationinfo)[1]
            djonesmatrices[1, 2, :, ant] = gentimeseries(obs.yamlconf["instrumentalpol"]["mode"], obs.stationinfo.d_pol1_loc[ant], obs.stationinfo.d_pol1_scale[ant], 0.0, obs.numchan, obs.rngcorrupt)
            djonesmatrices[2, 1, :, ant] = gentimeseries(obs.yamlconf["instrumentalpol"]["mode"], obs.stationinfo.d_pol2_loc[ant], obs.stationinfo.d_pol2_scale[ant], 0.0, obs.numchan, obs.rngcorrupt)

	    for chan in 1:obs.numchan
    	        if uppercase(obs.stationinfo.mount[ant]) == "ALT-AZ"
        	    polrotmatrices[1, 2, chan, :, ant] = djonesmatrices[1, 2, chan, ant] * exp.(2*im*(obs.stationinfo.feedangle_deg[ant].+parallacticanglematrix[:, ant]))
        	    polrotmatrices[2, 1, chan, :, ant] = djonesmatrices[2, 1, chan, ant] * exp.(-2*im*(obs.stationinfo.feedangle_deg[ant].+parallacticanglematrix[:, ant]))
		elseif uppercase(obs.stationinfo.mount[ant]) == "ALT-AZ+NASMYTH-L"
		    polrotmatrices[1, 2, chan, :, ant] = djonesmatrices[1, 2, chan, ant] * exp.(2*im*(obs.stationinfo.feedangle_deg[ant].+parallacticanglematrix[:, ant].-elevationmatrix[:, ant]))
        	    polrotmatrices[2, 1, chan, :, ant] = djonesmatrices[2, 1, chan, ant] * exp.(-2*im*(obs.stationinfo.feedangle_deg[ant].+parallacticanglematrix[:, ant].-elevationmatrix[:, ant]))
		elseif uppercase(obs.stationinfo.mount[ant]) == "ALT-AZ+NASMYTH-R"
		    polrotmatrices[1, 2, chan, :, ant] = djonesmatrices[1, 2, chan, ant] * exp.(2*im*(obs.stationinfo.feedangle_deg[ant].+parallacticanglematrix[:, ant].+elevationmatrix[:, ant]))
        	    polrotmatrices[2, 1, chan, :, ant] = djonesmatrices[2, 1, chan, ant] * exp.(-2*im*(obs.stationinfo.feedangle_deg[ant].+parallacticanglematrix[:, ant].+elevationmatrix[:, ant]))
		end
    	    end
        end

	# apply Dterms to data
	row = 1 # variable to index obs.data array
        for timeval in 1:size(uniqtimes)[1]
            for chan in 1:obs.numchan
                obs.data[:,:,chan,row] = polrotmatrices[:,:,chan,timeval,obs.antenna1[row]+1]*obs.data[:,:,chan,row]*adjoint(polrotmatrices[:,:,chan,timeval,obs.antenna2[row]+1])
            end
	    row += 1 # increment obs.data last index i.e. row number
        end
	
    else
        # compute in antenna frame
    end

    # write polarization matrices to HDF5 file
    g["djonesmatrices"] = djonesmatrices
    g["polrotmatrices"] = polrotmatrices
    attributes(g)["datatype"] = string(typeof(read(g[keys(g)[1]])))
    # close h5 file
    close(fid)

    @info("Compute and apply instrumental polarization... 🙆")
    
end
