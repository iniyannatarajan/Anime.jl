export stationgains

"""
    stationgains(h5file::String, scanno::Vector{Int}, times::Vector{Float64}, exposure::Float64, data::Array{Complex{Float32},4},
    stationinfo::DataFrame, mode::String, rngcorrupt::AbstractRNG, antenna1::Vector{Int}, antenna2::Vector{Int}, numchan::Int64)

Compute time-variable station gains and apply to data. The actual numerical values are serialized as HDF5.
"""
function stationgains(h5file::String, scanno::Vector{Int}, times::Vector{Float64}, exposure::Float64, data::Array{Complex{Float32},4},
    stationinfo::DataFrame, mode::String, rngcorrupt::AbstractRNG, antenna1::Vector{Int}, antenna2::Vector{Int}, numchan::Int64)
    # open h5 file for writing
    fid = h5open(h5file, "r+")
    g = create_group(fid, "stationgains")
    attributes(g)["desc"] = "Numerical values of time-variable per station G-Jones terms applied to data"
    attributes(g)["dims"] = "2 x 2 x ntimes_per_scan x nant" #for each scan, a 4d array of 2 x 2 x ntime x nant is stored

    # get unique scan numbers
    uniqscans = unique(scanno)

    # loop through each station to create a vector of 2x2 G-Jones terms evolving over time
    # TODO parametrise in terms of gain ratio
    row = 1 # variable to index data array
    for scan in uniqscans
        # compute ideal ntimes per scan
        actualtscanvec = unique(getindex(times, findall(scanno.==scan)))
	actualtscanveclen = length(actualtscanvec)
	idealtscanvec = collect(first(actualtscanvec):exposure:last(actualtscanvec))
	idealtscanveclen = length(idealtscanvec)

	# create 4d array to hold G-Jones terms per time per station
	gjonesmatrices = zeros(eltype(data), 2, 2, idealtscanveclen, size(stationinfo)[1]) # 2 x 2 x ntimes x nant

	for ant in eachindex(stationinfo.station)
        gjonesmatrices[1, 1, :, ant] = gentimeseries!(gjonesmatrices[1, 1, :, ant], mode, stationinfo.g_pol1_loc[ant], stationinfo.g_pol1_scale[ant], 0.0, idealtscanveclen, rngcorrupt)
	    gjonesmatrices[2, 2, :, ant] = gentimeseries!(gjonesmatrices[2, 2, :, ant], mode, stationinfo.g_pol2_loc[ant], stationinfo.g_pol2_scale[ant], 0.0, idealtscanveclen, rngcorrupt)
	end

        # loop over time/row and apply gjones terms corresponding to each baseline
	findnearest(A,x) = argmin(abs.(A .- x)) # define function to find nearest neighbour
        for t in 1:actualtscanveclen
            idealtimeindex = findnearest(idealtscanvec, actualtscanvec[t])
            # read all baselines present in a given time
            ant1vec = getindex(antenna1, findall(times.==actualtscanvec[t]))
            ant2vec = getindex(antenna2, findall(times.==actualtscanvec[t]))
            for (ant1,ant2) in zip(ant1vec, ant2vec)
                for chan in 1:numchan
                    data[:,:,chan,row] = gjonesmatrices[:,:,idealtimeindex,ant1+1]*data[:,:,chan,row]*adjoint(gjonesmatrices[:,:,idealtimeindex,ant2+1])
                end
                row += 1 # increment data last index i.e. row number
            end
        end

	# write to h5 file
	g["gjones_scan$(scan)"] = gjonesmatrices
    end

    # add datatype attribute
    attributes(g)["datatype"] = string(typeof(read(g[keys(g)[1]])))

    # close h5 file
    close(fid)

    @info("Compute and apply station gains 🙆")
end
