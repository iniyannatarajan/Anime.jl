export pointing

using Statistics

function compute_mispointvec(times::Vector{Float64}, pointinginterval::Float64, mispointsperscan::Int64)::Vector{Int32}
    """
    compute mispoint vector that holds mispoint id per scan for indexing into pointing amplitude error matrix
    """
    # derive some quantities
    ntimes = size(times)[1]
    mispointvec = zeros(Int32, ntimes)
    startval = first(times)
    stopval::Float64 = startval+pointinginterval
    mispointindex = 1

    # assign mispointvec values
    for ii in 1:ntimes
        if startval <= times[ii] <= stopval
            mispointvec[ii] = mispointindex
            continue
        end
        startval = stopval
        stopval = startval+pointinginterval
        mispointindex += 1
        mispointvec[ii] = mispointindex
    end

    return mispointvec
end

function longtermpointing()
end

function pointing(obs::CjlObservation)
    """
    Compute pointing errors and apply to data
    """
    @info("Computing pointing errors...")
    # get element type to be used
    elemtype = typeof(obs.data[1][1])
    nant = size(obs.stationinfo)[1]
 
    # open h5 file for writing
    fid = h5open(obs.yamlconf["hdf5corruptions"], "r+")
    g = create_group(fid, "pointingerrors")
    attributes(g)["desc"] = "Numerical values of time-variable short and long term pointing errors"
    attributes(g)["dims"] = "short and long term pointing errors have different dims"
 
    # get unique scan numbers
    uniqscans = unique(obs.scanno)

    pbfwhm = obs.stationinfo.pbfwhm230_arcsec./(mean(obs.chanfreqvec)/230.0e9) # scale primary beam to centre frequency of spw
    pointinginterval = obs.yamlconf["pointing"]["interval"] == "coherencetime" ? mean(obs.stationinfo.ctime_sec) : obs.yamlconf["pointing"]["interval"]
    @info("Generating new mispointings every $(pointinginterval) seconds")

    # TODO calculate antenna rise and set times and mask pointing offsets? Is this really necessary? 
    # If the source is not above horizon, the antenna would automatically be flagged; shouldn't make a difference
    row = 1 # variable to index obs.data array
    for scan in uniqscans
        # compute ideal ntimes per scan
        actualtscanvec = unique(getindex(obs.times, findall(obs.scanno.==scan)))
        actualtscanveclen = length(actualtscanvec)
        idealtscanvec = collect(first(actualtscanvec):obs.exposure:last(actualtscanvec))
        idealtscanveclen = length(idealtscanvec)

	mispointsperscan::Int64 = max(1, ceil((last(idealtscanvec)-first(idealtscanvec))/pointinginterval))
	perscanpointingoffsets = zeros(Float64, mispointsperscan, nant)
	perscanpointingamperrors = zeros(Float64, mispointsperscan, nant)

	# loop through stations and compute offsets and amplitude errors
	for ant in 1:nant
	    perscanpointingoffsets[:, ant] = gentimeseries("gaussian", 0.0f0, obs.stationinfo.pointingrms_arcsec[ant], 0.0, mispointsperscan, obs.rngcorrupt)
	    if obs.stationinfo.pbmodel[ant] == "gaussian"
                perscanpointingamperrors[:, ant] = exp.(-0.5.*(perscanpointingoffsets[:, ant]./(pbfwhm[ant]/2.35)).^2)
	    end
	end

	# compute mispointvec
	mispointvec = compute_mispointvec(idealtscanvec, pointinginterval, mispointsperscan)

	# TODO loop over data and apply pointing errors
        findnearest(A,x) = argmin(abs.(A .- x)) # define function to find nearest neighbour
        for t in 1:actualtscanveclen
            idealtimeindex = findnearest(idealtscanvec, actualtscanvec[t])
            # read all baselines present in a given time
            ant1vec = getindex(obs.antenna1, findall(obs.times.==actualtscanvec[t]))
            ant2vec = getindex(obs.antenna2, findall(obs.times.==actualtscanvec[t]))
            for (ant1,ant2) in zip(ant1vec, ant2vec)
                for chan in 1:obs.numchan
                    obs.data[:,:,chan,row] = perscanpointingamperrors[idealtimeindex,ant1+1]*obs.data[:,:,chan,row]*adjoint(perscanpointingamperrors[idealtimeindex,ant2+1])
                end
                row += 1 # increment obs.data last index i.e. row number
            end
        end

        # write to h5 file
	g["perscanamperr_scan$(scan)"] = perscanpointingamperrors
    end

    # add datatype attribute
    attributes(g)["datatype"] = string(typeof(read(g[keys(g)[1]])))

    # close h5 file
    close(fid)

    @info("Apply pointing errors to visibilities... 🙆")
end
