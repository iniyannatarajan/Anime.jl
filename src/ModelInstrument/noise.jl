export thermalnoise

"""
    thermalnoise(obs::CjlObservation)

Compute per-baseline thermal noise in visibility domain and apply to data. The actual numerical values are serialized as HDF5.
"""
function thermalnoise(obs::CjlObservation)
    # get type and size of the matrix to be created -- data is a 4d array
    elemtype = typeof(obs.data[1])

    # get unique times
    uniqtimes = unique(obs.times)
    ntimes = size(uniqtimes)[1]
    
    # open h5 file for writing
    fid = h5open(obs.yamlconf["hdf5corruptions"], "r+")
    g = create_group(fid, "thermalnoise")
    attributes(g)["desc"] = "Numerical values of thermal noise corruptions added to data"
    attributes(g)["dims"] = "stokes x nchan x ntimes_per_baseline (for all scans)"

    # get ant1 and ant2 vectors with unique elements
    uniqant1 = unique(obs.antenna1)
    uniqant2 = unique(obs.antenna2)

    # loop through each baseline pair and compute and apply thermal noise to obs.data
    thermalnoiserms = zeros(Float64, size(obs.data))
    thermalnoise = zeros(elemtype, size(obs.data))
    for a1 in uniqant1
	    for a2 in uniqant2
	        if a2>a1
  		        # compute sigma per baseline
                indices = intersect(findall(obs.antenna1.==a1), findall(obs.antenna2.==a2))
                thermalnoiserms[:, :, :, indices] .= sigmaperbl = (1/obs.yamlconf["correff"]) * sqrt((obs.stationinfo.sefd_Jy[a1+1]*obs.stationinfo.sefd_Jy[a2+1])
								/(2*obs.exposure*obs.chanwidth))
                # compute and add thermal noise to data
                thermalnoise[:, :, :, indices] = sigmaperbl*randn(obs.rngcorrupt, elemtype, size(obs.data[:,:,:,indices]))
	 	        obs.data[:, :, :, indices] = obs.data[:, :, :, indices] + thermalnoise[:, :, :, indices]
            end
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

    @info("Compute and apply thermal noise 🙆")
end
