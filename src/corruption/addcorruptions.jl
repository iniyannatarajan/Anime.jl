export addcorruptions

include(joinpath("troposphere.jl"))
#include(joinpath("ionosphere.jl"))
include(joinpath("beam.jl"))
include(joinpath("instrumentalpol.jl"))
include(joinpath("stationgains.jl"))
include(joinpath("bandpass.jl"))
include(joinpath("noise.jl"))

function computeweights!(totalrmsspec::Array{Float32, 4}, totalwtspec::Array{Float32, 4}, obs::CjlObservation)
    """
    Compute total rms (sigma) values from independent noise sources
    """
    fid = h5open(obs.yamlconf["hdf5corruptions"], "r")
    if haskey(fid, "thermalnoise") && haskey(fid["thermalnoise"], "thermalnoiserms")
        totalrmsspec[:,:,:,:] = read(fid["thermalnoise"]["thermalnoiserms"]).^2
        #@info(totalrmsspec[:,:,1,1], totalrmsspec[:,:,1,end])
    end
    if haskey(fid, "troposphere") && haskey(fid["troposphere"], "skynoiserms")
        totalrmsspec[:,:,:,:] = totalrmsspec[:,:,:,:] + read(fid["troposphere"]["skynoiserms"]).^2
    end

    # compute sigma_spectrum
    totalrmsspec[:,:,:,:] = sqrt.(totalrmsspec[:,:,:,:])
    #@info(totalrmsspec[:,:,1,1], totalrmsspec[:,:,1,end])

    # compute weight_spectrum
    totalwtspec[:,:,:,:] = 1 ./(totalrmsspec[:,:,:,:].^2)

    return totalrmsspec, totalwtspec
end

function addcorruptions(obs::CjlObservation)
    """
    Main fn to add corruptions
    """
    # create HDF5 file to store all corruptions
    fid = h5open(obs.yamlconf["hdf5corruptions"], "w") # using mode "w" to destroy existing contents
    close(fid)

    # add tropospheric effects
    obs.yamlconf["troposphere"]["enable"] && troposphere(obs)

    # add instrumental polarization
    obs.yamlconf["instrumentalpol"]["enable"] && @time instrumentalpol(obs)

    # add pointing errors
    obs.yamlconf["pointing"]["enable"] && @time pointing(obs)

    # add station gains
    obs.yamlconf["stationgains"]["enable"] && @time stationgains(obs)

    # add bandpasses
    obs.yamlconf["bandpass"]["enable"] && @time bandpass(obs)

    # add thermal noise
    obs.yamlconf["thermalnoise"]["enable"] && @time thermalnoise(obs)

    # replace NaNs with zeros
    obs.data[isnan.(obs.data)] .= 0.0+0.0*im

    # compute weight columns
    #@warn("WEIGHT_SPECTRUM and SIGMA_SPECTRUM are not filled in properly at the moment, while the code is being optimised.")
    totalrmsspec = zeros(Float32, 2, 2, obs.numchan, size(obs.data)[4]) # noise rms zero by default
    totalwtspec = ones(Float32, 2, 2, obs.numchan, size(obs.data)[4]) # weights unity by default

    @info("Computing weight and sigma spectrum arrays...")
    @time computeweights!(totalrmsspec, totalwtspec, obs)

    # convert the sigma_spec and weight_spec arrays to the format required by Casacore.jl
    revtotalrmsspec3dres = permutedims(totalrmsspec, (2,1,3,4))
    revtotalrmsspec3d = reshape(revtotalrmsspec3dres, 4, size(revtotalrmsspec3dres)[3], :)
    revtotalrmsspec = [Matrix{Float32}(revtotalrmsspec3d[:,:,i]) for i in 1:size(revtotalrmsspec3d)[3]]

    revtotalwtspec3dres = permutedims(totalwtspec, (2,1,3,4))
    revtotalwtspec3d = reshape(revtotalwtspec3dres, 4, size(revtotalwtspec3dres)[3], :)
    revtotalwtspec = [Matrix{Float32}(revtotalwtspec3d[:,:,i]) for i in 1:size(revtotalwtspec3d)[3]]

    # compute channel averaged weight and sigma column values
    revtotalrms = Vector{Vector{Float32}}()
    revtotalwt = Vector{Vector{Float32}}()
    for row in 1:size(obs.data)[4] # nrows
        tmpvec = dropdims(mean(revtotalrmsspec[row], dims=2), dims=2)
        push!(revtotalrms, tmpvec) # average the sigma values
        push!(revtotalwt, 1 ./tmpvec.^2)
    end

    # reshape data array for Casacore.jl
    revdata3dres = permutedims(obs.data, (2,1,3,4))
    revdata3d = reshape(revdata3dres, 4, size(revdata3dres)[3], :)
    revdata = [Matrix{ComplexF32}(revdata3d[:,:,i]) for i in 1:size(revdata3d)[3]]

    # when all the corruptions have been applied, write the above columns to ms
    table = CCTable(obs.yamlconf["msname"], CCTables.Update)
    table[:DATA] = revdata # Float32 to conform to the MSv2 specification (which WSClean expects... sometimes!)
    table[:SIGMA_SPECTRUM] = revtotalrmsspec
	table[:WEIGHT_SPECTRUM] = revtotalwtspec
    table[:SIGMA] = revtotalrms
	table[:WEIGHT] = revtotalwt

    @info("Write arrays to disk... 🙆")
end
