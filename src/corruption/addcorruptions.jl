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

    # add thermal noise
    obs.yamlconf["thermalnoise"]["enable"] && thermalnoise(obs)
    
    # when all the corruptions have been applied, write data column back to ms
    table = CCTable(obs.yamlconf["msname"], CCTables.Update)
    table[:DATA] = obs.data # Float32 to conform to the MSv2 specification (which WSClean expects... sometimes!)
    @info("Corrupted data written to Measurement Set 🙆")
end
