export addcorruptions

include(joinpath("troposphere.jl"))
include(joinpath("ionosphere.jl"))
include(joinpath("beam.jl"))
include(joinpath("instrumentpol.jl"))
include(joinpath("stationgains.jl"))
include(joinpath("bandpass.jl"))
include(joinpath("noise.jl"))

function addcorruptions(obs::CjlObservation)
    """
    Main fn to add corruptions
    """
    # create HDF5 file to store all corruptions
    fid = h5open(obs.yamlconf["corrupth5name"], "w")
    close(fid)

    # add station gains
    obs.yamlconf["stationgains"]["enable"] && stationgains(obs)

    # add thermal noise
    obs.yamlconf["thermalnoise"]["enable"] && thermalnoise(obs)
   
    # convert data array back to the format required by Casacore
    revdata3dres = permutedims(obs.data, (2,1,3,4))
    revdata3d = reshape(revdata3dres, 4, size(revdata3dres)[3], :)
    revdata = [Matrix{ComplexF32}(revdata3d[:,:,i]) for i in 1:size(revdata3d)[3]]

    # when all the corruptions have been applied, write data column back to ms
    table = CCTable(obs.yamlconf["msname"], CCTables.Update)
    table[:DATA] = revdata # Float32 to conform to the MSv2 specification (which WSClean expects... sometimes!)
    @info("Corrupted data written to Measurement Set 🙆")
end
