export instrumentalpol

"""
    instrumentalpol(obs::CjlObservation; h5file::String="")

Compute instrumental polarization (leakage, or "D-Jones" terms) and apply to data. The actual numerical values are serialized as HDF5.
"""
function instrumentalpol(obs::CjlObservation; h5file::String="")
    # get unique scan numbers
    uniqscans = unique(obs.scanno)

    # get unique times
    uniqtimes = unique(obs.times)
    ntimes = size(uniqtimes)[1]

    # compute necessary quantities
    elevationmatrix = elevationangle(obs.times, obs.phasedir, obs.stationinfo, obs.pos)
    parallacticanglematrix = parallacticangle(obs)

	# D-terms -- perform twice the feed angle rotation
	djonesmatrices = ones(eltype(obs.data), 2, 2, obs.numchan, size(obs.stationinfo)[1]) # 2 x 2 x nchan x nant
	polrotmatrices = ones(eltype(obs.data), 2, 2, obs.numchan, ntimes, size(obs.stationinfo)[1])

    # compute D-terms
    if obs.yamlconf["instrumentalpol"]["visibilityframe"] == "antenna"
        # compute in antenna frame
        @warn("Visibilities in antenna frame not implemented yet. Skipping instrumental polarization... ")
    else
        if obs.yamlconf["instrumentalpol"]["visibilityframe"] != "sky"
	        @warn("Visibility frame set to $(obs.yamlconf["instrumentalpol"]["visibilityframe"]). Computing visibilities in 'sky' frame...")
        else
            @info("Applying instrumental polarization and recording visibilities in 'sky' frame...")
        end

	    for ant in eachindex(obs.stationinfo.station)
            djonesmatrices[1, 2, :, ant] = gentimeseries!(djonesmatrices[1, 2, :, ant], obs.yamlconf["instrumentalpol"]["mode"], obs.stationinfo.d_pol1_loc[ant], obs.stationinfo.d_pol1_scale[ant], 0.0, obs.numchan, obs.rngcorrupt)
            djonesmatrices[2, 1, :, ant] = gentimeseries!(djonesmatrices[2, 1, :, ant], obs.yamlconf["instrumentalpol"]["mode"], obs.stationinfo.d_pol2_loc[ant], obs.stationinfo.d_pol2_scale[ant], 0.0, obs.numchan, obs.rngcorrupt)

	        for t in 1:ntimes
    	        for chan in 1:obs.numchan
    	            if uppercase(obs.stationinfo.mount[ant]) == "ALT-AZ"
		                polrotmatrices[1, 2, chan, t, ant] = djonesmatrices[1, 2, chan, ant] * exp(2*im*(deg2rad(obs.stationinfo.feedangle_deg[ant])+parallacticanglematrix[t, ant]))
		                polrotmatrices[2, 1, chan, t, ant] = djonesmatrices[2, 1, chan, ant] * exp(-2*im*(deg2rad(obs.stationinfo.feedangle_deg[ant])+parallacticanglematrix[t, ant]))
		            elseif uppercase(obs.stationinfo.mount[ant]) == "ALT-AZ+NASMYTH-L"
    		            polrotmatrices[1, 2, chan, t, ant] = djonesmatrices[1, 2, chan, ant] * exp(2*im*(deg2rad(obs.stationinfo.feedangle_deg[ant])+parallacticanglematrix[t, ant]-elevationmatrix[t, ant]))
		                polrotmatrices[2, 1, chan, t, ant] = djonesmatrices[2, 1, chan, ant] * exp(-2*im*(deg2rad(obs.stationinfo.feedangle_deg[ant])+parallacticanglematrix[t, ant]-elevationmatrix[t, ant]))
		            elseif uppercase(obs.stationinfo.mount[ant]) == "ALT-AZ+NASMYTH-R"
		                polrotmatrices[1, 2, chan, t, ant] = djonesmatrices[1, 2, chan, ant] * exp(2*im*(deg2rad(obs.stationinfo.feedangle_deg[ant])+parallacticanglematrix[t, ant]+elevationmatrix[t, ant]))
		                polrotmatrices[2, 1, chan, t, ant] = djonesmatrices[2, 1, chan, ant] * exp(-2*im*(deg2rad(obs.stationinfo.feedangle_deg[ant])+parallacticanglematrix[t, ant]+elevationmatrix[t, ant]))
		            end
    	        end
	        end
        end

	    # apply Dterms to data TODO rewrite the following loop
        row = 1 # variable to index obs.data array
        for scan in uniqscans
            # compute ideal ntimes per scan
            actualtscanvec = unique(getindex(obs.times, findall(obs.scanno.==scan)))
            actualtscanveclen = length(actualtscanvec)
            idealtscanvec = collect(first(actualtscanvec):obs.exposure:last(actualtscanvec))
            idealtscanveclen = length(idealtscanvec)
    
            # loop over time/row and apply Dterms corresponding to each baseline
            findnearest(A,x) = argmin(abs.(A .- x)) # define function to find nearest neighbour
            for t in 1:actualtscanveclen
                idealtimeindex = findnearest(idealtscanvec, actualtscanvec[t])
                # read all baselines present in a given time
                ant1vec = getindex(obs.antenna1, findall(obs.times.==actualtscanvec[t]))
                ant2vec = getindex(obs.antenna2, findall(obs.times.==actualtscanvec[t]))
                for (ant1,ant2) in zip(ant1vec, ant2vec)
                    for chan in 1:obs.numchan
                        obs.data[:,:,chan,row] = polrotmatrices[:,:,chan,idealtimeindex,ant1+1]*obs.data[:,:,chan,row]*adjoint(polrotmatrices[:,:,chan,idealtimeindex,ant2+1])
                    end
                    row += 1 # increment obs.data last index i.e. row number
                end
            end
        end
    end

    # write polarization matrices to HDF5 file
    if !isempty(h5file)
        fid = h5open(h5file, "cw")
        g = create_group(fid, "polarization")
        attributes(g)["desc"] = "Numerical values of instrumental polarization matrices applied to data (Dterms and fullpolproduct)"
        attributes(g)["dims"] = "2 x 2 x nchannels x nant" #for each scan, a 4d array of 2 x 2 x nchan x nant is stored

        if !(haskey(g, "elevation"))
            g["elevation"] = elevationmatrix
        end

        if !(haskey(g, "parallacticangle"))
            g["parallacticangle"] = parallacticanglematrix
        end

        g["djonesmatrices"] = djonesmatrices
        g["polrotmatrices"] = polrotmatrices
        attributes(g)["datatype"] = string(typeof(read(g[keys(g)[1]])))
        close(fid)
    end

    @info("Compute and apply instrumental polarization... 🙆")
    
end
