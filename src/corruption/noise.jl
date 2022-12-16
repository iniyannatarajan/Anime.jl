export thermalnoise

function thermalnoise(obs::CjlObservation)
    """
    Add thermal noise to visibilities
    """
    # get matrix type and size to be created -- data is a 4d array
    elemtype = typeof(obs.data[1])

    # get unique times
    uniqtimes = unique(obs.times)
    ntimes = size(uniqtimes)[1]
    
    # open h5 file for writing
    fid = h5open(obs.yamlconf["hdf5corruptions"], "r+")
    g = create_group(fid, "thermalnoise")
    attributes(g)["desc"] = "Numerical values of thermal noise corruptions added to data"
    attributes(g)["dims"] = "stokes x nchan x ntimes_per_baseline (for all scans)"

    #=
    # get ant1 and ant2 vectors with unique elements
    uniqant1 = unique(obs.antenna1)
    uniqant2 = unique(obs.antenna2)

    # loop through each baseline pair and compute and apply thermal noise to obs.data
    for a1 in uniqant1
	for a2 in uniqant2
	    if a2>a1
		# compute sigma per baseline
		sigmaperbl = (1/obs.yamlconf["correff"]) * sqrt((obs.stationinfo.sefd_Jy[a1+1]*obs.stationinfo.sefd_Jy[a2+1])
								/(2*obs.exposure*obs.chanwidth))
		indices = intersect(findall(obs.antenna1.==a1), findall(obs.antenna2.==a2))

                # compute and add thermal noise to data
		thermalvec = sigmaperbl*randn(obs.rngcorrupt, elemtype, size(obs.data[:,:,:,indices]))
		obs.data[:,:,:,indices] += thermalvec

                # write as individual dataset within the group created above in the h5 file
                g["baseline_$(obs.stationinfo.station[a1+1])-$(obs.stationinfo.station[a2+1])"] = thermalvec #reduce((x,y) -> cat(x, y, dims=3), thermalvec)
	    end
	end
    end
    =#

    # compute thermal noise and add to data
    thermalnoiserms = zeros(Float64, size(obs.data))
    thermalnoise = zeros(elemtype, size(obs.data))
    row = 1
    for t in 1:ntimes # no. of unique times
        # read all baselines present in a given time
        ant1vec = getindex(obs.antenna1, findall(obs.times.==uniqtimes[t]))
        ant2vec = getindex(obs.antenna2, findall(obs.times.==uniqtimes[t]))
        for (ant1,ant2) in zip(ant1vec, ant2vec)
            for chan in 1:obs.numchan
		thermalnoiserms[:, :, chan, row] .= (1/obs.yamlconf["correff"]) * sqrt((obs.stationinfo.sefd_Jy[ant1+1]*obs.stationinfo.sefd_Jy[ant2+1])/(2*obs.exposure*obs.chanwidth))
		for jj in 1:2
		    for ii in 1:2
                        thermalnoise[ii, jj, chan, row] = thermalnoiserms[ii, jj, chan, row]*randn(obs.rngcorrupt, Float64) # thermal noise is polarized
                        obs.data[ii, jj, chan, row] += thermalnoise[ii, jj, chan, row]
		    end
	        end
            end
            row += 1 # increment the last dimension i.e. row number
        end
    end

    # write to h5 file
    if !(haskey(g, "thermalnoiserms"))
        g["thermalnoiserms"] = thermalnoiserms
    end
    if !(haskey(g, "thermalnoise"))
        g["thermalnoise"] = thermalnoise
    end

    # add datatype attribute
    attributes(g)["datatype"] = string(typeof(read(g[keys(g)[1]])))

    # close h5 file
    close(fid)

    @info("Compute and apply thermal noise... 🙆")
end
