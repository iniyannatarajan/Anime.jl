export instrumentalpol

"""
    instrumentalpol(scanno::Vector{Int32}, times::Vector{Float64}, stationinfo::DataFrame, phasedir::Array{Float64,2},
    pos::Array{Float64, 2}, data::Array{Complex{Float32},4}, numchan::Int64, polframe::String, polmode::String,
    antenna1::Vector{Int32}, antenna2::Vector{Int32}, exposure::Float64, rngcorrupt::AbstractRNG; h5file::String="")

Compute instrumental polarization (leakage, or "D-Jones" terms) and apply to data. The actual numerical values are serialized as HDF5.
"""
function instrumentalpol(scanno::Vector{Int32}, times::Vector{Float64}, stationinfo::DataFrame, phasedir::Array{Float64,2},
    pos::Array{Float64, 2}, data::Array{Complex{Float32},4}, numchan::Int64, polframe::String, polmode::String,
    antenna1::Vector{Int32}, antenna2::Vector{Int32}, exposure::Float64, rngcorrupt::AbstractRNG; h5file::String="")
    # get unique scan numbers
    uniqscans = unique(scanno)

    # get unique times
    uniqtimes = unique(times)
    ntimes = size(uniqtimes)[1]

    # compute necessary quantities
    elevationmatrix = elevationangle(times, phasedir, stationinfo, pos)
    parallacticanglematrix = parallacticangle(times, phasedir, stationinfo, pos)

	# D-terms -- perform twice the feed angle rotation
	djonesmatrices = ones(eltype(data), 2, 2, numchan, size(stationinfo)[1]) # 2 x 2 x nchan x nant
	polrotmatrices = ones(eltype(data), 2, 2, numchan, ntimes, size(stationinfo)[1])

    # compute D-terms
    if polframe == "antenna"
        # compute in antenna frame
        @warn("Visibilities in antenna frame not implemented yet. Skipping instrumental polarization... ")
    else
        if polframe != "sky"
	        @warn("Visibility frame set to $(polframe). Computing visibilities in 'sky' frame...")
        else
            @info("Applying instrumental polarization and recording visibilities in 'sky' frame...")
        end

	    for ant in eachindex(stationinfo.station)
            djonesmatrices[1, 2, :, ant] = gentimeseries!(djonesmatrices[1, 2, :, ant], polmode, stationinfo.d_pol1_loc[ant], Float32(stationinfo.d_pol1_scale[ant]), Float32(0.0), numchan, rngcorrupt)
            djonesmatrices[2, 1, :, ant] = gentimeseries!(djonesmatrices[2, 1, :, ant], polmode, stationinfo.d_pol2_loc[ant], Float32(stationinfo.d_pol2_scale[ant]), Float32(0.0), numchan, rngcorrupt)

	        for t in 1:ntimes
    	        for chan in 1:numchan
    	            if uppercase(stationinfo.mount[ant]) == "ALT-AZ"
		                polrotmatrices[1, 2, chan, t, ant] = djonesmatrices[1, 2, chan, ant] * exp(2*im*(deg2rad(stationinfo.feedangle_deg[ant])+parallacticanglematrix[t, ant]))
		                polrotmatrices[2, 1, chan, t, ant] = djonesmatrices[2, 1, chan, ant] * exp(-2*im*(deg2rad(stationinfo.feedangle_deg[ant])+parallacticanglematrix[t, ant]))
		            elseif uppercase(stationinfo.mount[ant]) == "ALT-AZ+NASMYTH-L"
    		            polrotmatrices[1, 2, chan, t, ant] = djonesmatrices[1, 2, chan, ant] * exp(2*im*(deg2rad(stationinfo.feedangle_deg[ant])+parallacticanglematrix[t, ant]-elevationmatrix[t, ant]))
		                polrotmatrices[2, 1, chan, t, ant] = djonesmatrices[2, 1, chan, ant] * exp(-2*im*(deg2rad(stationinfo.feedangle_deg[ant])+parallacticanglematrix[t, ant]-elevationmatrix[t, ant]))
		            elseif uppercase(stationinfo.mount[ant]) == "ALT-AZ+NASMYTH-R"
		                polrotmatrices[1, 2, chan, t, ant] = djonesmatrices[1, 2, chan, ant] * exp(2*im*(deg2rad(stationinfo.feedangle_deg[ant])+parallacticanglematrix[t, ant]+elevationmatrix[t, ant]))
		                polrotmatrices[2, 1, chan, t, ant] = djonesmatrices[2, 1, chan, ant] * exp(-2*im*(deg2rad(stationinfo.feedangle_deg[ant])+parallacticanglematrix[t, ant]+elevationmatrix[t, ant]))
		            end
    	        end
	        end
        end

	    # apply Dterms to data TODO rewrite the following loop
        row = 1 # variable to index data array
        for scan in uniqscans
            # compute ideal ntimes per scan
            actualtscanvec = unique(getindex(times, findall(scanno.==scan)))
            actualtscanveclen = length(actualtscanvec)
            idealtscanvec = collect(first(actualtscanvec):exposure:last(actualtscanvec))
            idealtscanveclen = length(idealtscanvec)
    
            # loop over time/row and apply Dterms corresponding to each baseline
            findnearest(A,x) = argmin(abs.(A .- x)) # define function to find nearest neighbour
            for t in 1:actualtscanveclen
                idealtimeindex = findnearest(idealtscanvec, actualtscanvec[t])
                # read all baselines present in a given time
                ant1vec = getindex(antenna1, findall(times.==actualtscanvec[t]))
                ant2vec = getindex(antenna2, findall(times.==actualtscanvec[t]))
                for (ant1,ant2) in zip(ant1vec, ant2vec)
                    for chan in 1:numchan
                        data[:,:,chan,row] = polrotmatrices[:,:,chan,idealtimeindex,ant1+1]*data[:,:,chan,row]*adjoint(polrotmatrices[:,:,chan,idealtimeindex,ant2+1])
                    end
                    row += 1 # increment data last index i.e. row number
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
