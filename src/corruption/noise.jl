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
    for a1 in uniqant1
        for a2 in uniqant2
	    if a2>a1
		sigmaperbl = (1/obs.yamlconf["correff"]) * sqrt((obs.stationinfo.sefd_Jy[a1+1]*obs.stationinfo.sefd_Jy[a2+1])
								/(2*obs.yamlconf["inttime"]*(obs.yamlconf["spw"]["bandwidth"][1]*1e9/obs.yamlconf["spw"]["channels"][1])))
		indices = intersect(findall(obs.antenna1.==a1), findall(obs.antenna2.==a2))
                datasubset = getindex(obs.data, intersect(findall(obs.antenna1.==a1), findall(obs.antenna2.==a2)))

		# create thermalnoise vector of matrices
		thermalvec = mattype[]
		for ii in 1:length(datasubset)
		    push!(thermalvec, sigmaperbl*randn(obs.rngcorrupt, elemtype, matsize))
		end
                setindex!(obs.data, datasubset+thermalvec, indices)
	    end
	end
    end
    @info("Thermal noise applied 🙆")
end
