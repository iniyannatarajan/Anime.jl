export thermalnoise

function thermalnoise(obs::CjlObservation)
    """
    Add thermal noise to visibilities
    """
    # get matrix type and size to be created -- data is a 4d array
    elemtype = typeof(obs.data[1])

    # open h5 file for writing
    fid = h5open(obs.yamlconf["hdf5corruptions"], "r+")
    g = create_group(fid, "thermalnoise")
    attributes(g)["desc"] = "Numerical values of thermal noise corruptions added to data"
    attributes(g)["format"] = "each dataset corresponds to the 3d array stokes X channel X rows_per_baseline (for all scans)"

    # get ant1 and ant2 vectors with unique elements
    uniqant1 = unique(obs.antenna1)
    uniqant2 = unique(obs.antenna2)

    # loop through each baseline pair and compute and apply thermal noise to obs.data
    for a1 in uniqant1
	for a2 in uniqant2
	    if a2>a1
		# compute sigma per baseline
		sigmaperbl = (1/obs.yamlconf["thermalnoise"]["correff"]) * sqrt((obs.stationinfo.sefd_Jy[a1+1]*obs.stationinfo.sefd_Jy[a2+1])
								/(2*obs.exposure*obs.chanwidth))
		indices = intersect(findall(obs.antenna1.==a1), findall(obs.antenna2.==a2))

                # compute and add thermal noise to data
		thermalvec = sigmaperbl*randn(obs.rngcorrupt, elemtype, size(obs.data[:,:,:,indices]))
		obs.data[:,:,:,indices] += thermalvec

                # write as individual dataset within the group created above in the h5 file
                g["baseline $(obs.stationinfo.station[a1+1])-$(obs.stationinfo.station[a2+1])"] = thermalvec #reduce((x,y) -> cat(x, y, dims=3), thermalvec)
	    end
	end
    end

    # add datatype attribute
    attributes(g)["datatype"] = string(typeof(read(g[keys(g)[1]])))

    # close h5 file
    close(fid)

    @info("Compute and apply thermal noise... 🙆")
end
