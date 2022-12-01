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
    #= # get matrix type and size to be created -- data is a 4d array
    matsize = size(obs.data[1])
    mattype = typeof(obs.data[1])
    elemtype = typeof(obs.data[1][1])

    # open h5 file for writing
    fid = h5open(obs.yamlconf["hd5corruptions"], "r+")
    g = create_group(fid, "stationgains")
    attributes(g)["desc"] = "Numerical values of time-variable per station G-Jones terms applied to data"
    attributes(g)["format"] = "each dataset corresponds to the 5d array nscans x nant x ntimes x 2 x 2"=#

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
	gjonesmatrices = zeros(ComplexF32, 2, 2, idealtscanveclen, size(obs.stationinfo)[1]) # 2 x 2 x ntimes x nant

	for ant in 1:size(obs.stationinfo)[1]
	   gjonesmatrices[1, 1, :, ant] = generatetimeseries(obs.yamlconf["stationgains"]["mode"], obs.stationinfo.g_pol1_loc[ant], obs.stationinfo.g_pol1_scale[ant], idealtscanveclen)
	   gjonesmatrices[2, 2, :, ant] = generatetimeseries(obs.yamlconf["stationgains"]["mode"], obs.stationinfo.g_pol2_loc[ant], obs.stationinfo.g_pol2_scale[ant], idealtscanveclen)
	end

	# add to gjonesdict
	gjonesdict[scan] = gjonesmatrices

        # loop over time/row and apply gjones terms corresponding to each baseline
	for idealtimeindex in 1:idealtscanveclen
	    if idealtscanvec[idealtimeindex] in actualtscanvec
		for chan in 1:obs.numchan
		    obs.data[:,:,chan,idealtimeindex] = gjonesdict[scan][:,:,idealtimeindex,obs.antenna1[idealtimeindex]+1]*obs.data[:,:,chan,idealtimeindex]*adjoint(gjonesdict[scan][:,:,idealtimeindex,obs.antenna2[idealtimeindex]+1])
	        end
	    else
		continue
	    end
        end
    end

    #= # get ant1 and ant2 vectors with unique elements
    uniqant1 = unique(obs.antenna1)
    uniqant2 = unique(obs.antenna2)

    # loop through each baseline pair and compute and apply thermal noise to obs.data
    for a1 in uniqant1
	for a2 in uniqant2
	    if a2>a1
		# compute sigma per baseline
		sigmaperbl = (1/obs.yamlconf["correff"]) * sqrt((obs.stationinfo.sefd_Jy[a1+1]*obs.stationinfo.sefd_Jy[a2+1])
								/(2*obs.yamlconf["inttime"]*(obs.yamlconf["spw"]["bandwidth"][1]*1e9/obs.yamlconf["spw"]["channels"][1])))
		indices = intersect(findall(obs.antenna1.==a1), findall(obs.antenna2.==a2))
                datasubset = getindex(obs.data, indices)
                thermalvec = mattype[] # create empty thermalvec

		# write to thermalnoise vector of matrices
		for ii in 1:length(datasubset)
		    push!(thermalvec, sigmaperbl*randn(obs.rngcorrupt, elemtype, matsize))
		end

                # add thermal noise to data
                setindex!(obs.data, datasubset+thermalvec, indices)

                # reduce vector of matrices to 3d array for HDF5 compatibility and write as dataset within the group in the h5 file
                g["baseline $(a1)-$(a2)"] = reduce((x,y) -> cat(x, y, dims=3), thermalvec)
	    end
	end
    end

    # add datatype attribute
    attributes(g)["datatype"] = string(typeof(read(g[keys(g)[1]])))

    # close h5 file
    close(fid)=#

    @info("Compute and apply station gains... 🙆")
end
