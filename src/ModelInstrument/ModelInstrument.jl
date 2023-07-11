module ModelInstrument

using CSV
using HDF5
using Distributions
using DataFrames
using Casacore.Tables: Tables as CCTables, Table as CCTable

using ..Observe: CjlObservation
using ..Utils: parallacticangle, elevationangle, gentimeseries!

include(joinpath(@__DIR__, "troposphere.jl"))
include(joinpath(@__DIR__, "ionosphere.jl"))
include(joinpath(@__DIR__, "beam.jl"))
include(joinpath(@__DIR__, "instrumentalpol.jl"))
include(joinpath(@__DIR__, "stationgains.jl"))
include(joinpath(@__DIR__, "bandpass.jl"))
include(joinpath(@__DIR__, "noise.jl"))
include(joinpath(@__DIR__, "applymodel.jl"))

end