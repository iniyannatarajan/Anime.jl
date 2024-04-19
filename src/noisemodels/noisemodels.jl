export thermalnoise!

"""
    thermalnoise!(data::Array{Complex{Float32},4}, times::Vector{Float64}, antenna1::Vector{Int32}, antenna2::Vector{Int32},
    correff::Float64, exposure::Float64, chanwidth::Float64, corruptseed::Int, sefd::Vector{Float64}; h5file::String="",
    noisefile::String="")

Compute per-baseline thermal noise using radiometer equation and apply to data. The actual numerical values are serialized in HDF5 format.
"""
function thermalnoise!(data::Array{Complex{Float32},4}, times::Vector{Float64}, antenna1::Vector{Int32}, antenna2::Vector{Int32},
    correff::Float64, exposure::Float64, chanwidth::Float64, corruptseed::Int, sefd::Vector{Float64}; h5file::String="",
    noisefile::String="")
    @info("Computing thermal noise...")
    #initialize RNG with seed
    rngcorrupt = Xoshiro(corruptseed)

    # get unique times
    uniqtimes = unique(times)
    ntimes = size(uniqtimes)[1]

    # get ant1 and ant2 vectors with unique elements
    uniqant1 = unique(antenna1)
    uniqant2 = unique(antenna2)

    # loop through each baseline pair and compute and apply thermal noise to data
    thermalnoiserms = zeros(Float64, size(data))
    thermalnoise = zeros(eltype(data), size(data))

    isnoisetablevalid = false
    if isfile(noisefile)
        fidnoise = h5open(noisefile, "r")
        if size(read(fidnoise["thermalnoise"]["thermalnoiserms"])) != size(thermalnoiserms) || size(read(fidnoise["thermalnoise"]["thermalnoise"])) != size(thermalnoise) || eltype(read(fidnoise["thermalnoise"]["thermalnoiserms"])) != eltype(thermalnoiserms) || eltype(read(fidnoise["thermalnoise"]["thermalnoise"])) != eltype(thermalnoise)
            @warn("Dims/type mismatch between input $(noisefile) and observation parameters. Will compute thermal noise...")
        else
            thermalnoiserms = deepcopy(read(fidnoise["thermalnoise"]["thermalnoiserms"]))
            thermalnoise = deepcopy(read(fidnoise["thermalnoise"]["thermalnoise"]))
            isnoisetablevalid = true
        end
        close(fidnoise)
    end

    if isempty(noisefile) || !isnoisetablevalid 
        for a1 in uniqant1
            for a2 in uniqant2
                if a2>a1
 	                # compute sigma per baseline
                    indices = intersect(findall(antenna1.==a1), findall(antenna2.==a2))
                    thermalnoiserms[:, :, :, indices] .= sigmaperbl = (1/correff) * sqrt((sefd[a1+1]*sefd[a2+1])
	    			    			/(2*exposure*chanwidth))
                    # compute and add thermal noise to data
                    thermalnoise[:, :, :, indices] = sigmaperbl*randn(rngcorrupt, eltype(data), size(data[:,:,:,indices]))
 	 	            data[:, :, :, indices] = data[:, :, :, indices] + thermalnoise[:, :, :, indices]
                end
            end
        end
    end

    # write to h5 file
    if !isempty(h5file)
        fid = h5open(h5file, "cw")
        g = create_group(fid, "thermalnoise")
        HDF5.attributes(g)["desc"] = "Numerical values of thermal noise corruptions added to data"
        HDF5.attributes(g)["dims"] = "stokes x nchan x ntimes_per_baseline (for all scans)"

        if !(haskey(g, "thermalnoiserms"))
            g["thermalnoiserms"] = thermalnoiserms
        end
        if !(haskey(g, "thermalnoise"))
            g["thermalnoise"] = thermalnoise
        end

        # add datatype attribute
        HDF5.attributes(g)["datatype"] = string(typeof(read(g[keys(g)[1]])))

        # close h5 file
        close(fid)
    end

    @info("Compute and apply thermal noise 🙆")
end

"""
    thermalnoise!(ms::MeasurementSet, stationinfo::DataFrame, obsconfig::Dict; h5file::String="", noisefile::String="")

Alias for use in pipelines.
"""
function thermalnoise!(ms::MeasurementSet, stationinfo::DataFrame, obsconfig::Dict; h5file::String="", noisefile::String="")
    thermalnoise!(ms.data, ms.times, ms.antenna1, ms.antenna2, obsconfig["correff"], ms.exposure, ms.chanwidth,
     obsconfig["corruptseed"], stationinfo.sefd_Jy, h5file=h5file, noisefile=noisefile)
end