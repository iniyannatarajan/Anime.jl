export thermalnoise

"""
    thermalnoise(times::Vector{Float64}, antenna1::Vector{Int}, antenna2::Vector{Int}, data::Array{Complex{Float32},4}, correff::Float64,
    exposure::Float64, chanwidth::Float64, rngcorrupt::AbstractRNG, sefd::Vector{Float64}; h5file::String="")

Compute per-baseline thermal noise in visibility domain and apply to data. The actual numerical values are serialized as HDF5.
"""
function thermalnoise(times::Vector{Float64}, antenna1::Vector{Int}, antenna2::Vector{Int}, data::Array{Complex{Float32},4}, correff::Float64,
    exposure::Float64, chanwidth::Float64, rngcorrupt::AbstractRNG, sefd::Vector{Float64}; h5file::String="")
    # get unique times
    uniqtimes = unique(times)
    ntimes = size(uniqtimes)[1]

    # get ant1 and ant2 vectors with unique elements
    uniqant1 = unique(antenna1)
    uniqant2 = unique(antenna2)

    # loop through each baseline pair and compute and apply thermal noise to data
    thermalnoiserms = zeros(Float64, size(data))
    thermalnoise = zeros(eltype(data), size(data))
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

    # write to h5 file
    if !isempty(h5file)
        fid = h5open(h5file, "cw")
        g = create_group(fid, "thermalnoise")
        attributes(g)["desc"] = "Numerical values of thermal noise corruptions added to data"
        attributes(g)["dims"] = "stokes x nchan x ntimes_per_baseline (for all scans)"

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
    end

    @info("Compute and apply thermal noise 🙆")
end
