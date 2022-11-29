export thermalnoise

function thermalnoise(obs::CjlObservation)
    """
    Add thermal noise to visibilities
    """
    # get matrix type and size to be created
    matsize = size(obs.data[1])
    mattype = typeof(obs.data[1])
    elemtype = typeof(obs.data[1][1])

    uniqant1 = unique(obs.antenna1)
    uniqant2 = unique(obs.antenna2)
    thermalvec = mattype[] # create empty thermalvec
    for a1 in uniqant1
	for a2 in uniqant2
	    if a2>a1
		sigmaperbl = (1/obs.yamlconf["correff"]) * sqrt((obs.stationinfo.sefd_Jy[a1+1]*obs.stationinfo.sefd_Jy[a2+1])
								/(2*obs.yamlconf["inttime"]*(obs.yamlconf["spw"]["bandwidth"][1]*1e9/obs.yamlconf["spw"]["channels"][1])))
		indices = intersect(findall(obs.antenna1.==a1), findall(obs.antenna2.==a2))
                datasubset = getindex(obs.data, intersect(findall(obs.antenna1.==a1), findall(obs.antenna2.==a2)))

		# write to thermalnoise vector of matrices
		for ii in 1:length(datasubset)
		    push!(thermalvec, sigmaperbl*randn(obs.rngcorrupt, elemtype, matsize))
		end


	    end
	end
    end

    # add thermal noise to data
    setindex!(obs.data, obs.data+thermalvec, 1:length(thermalvec))

    # save thermal noise vector
    # reduce vector of matrices to 3d array for HDF5 compatibility
    thermalnoisevector = reduce((x,y) -> cat(x, y, dims=3), thermalvec)

    h5open(obs.yamlconf["corrupth5name"], "r+") do fid
        #=if haskey(fid, "thermalnoise")
       	    delete_object(fid, "thermalnoise")
        end=#
        g = create_group(fid, "thermalnoise")
        g["thermalnoisevector"] = thermalnoisevector
    end
    @info("Thermal noise applied 🙆")
end
