using Anime
using Test
using YAML
using Random

@testset "Load Observation" begin
    include("testloadobs.jl")
end

@testset "Instrument Models" begin
    include("testinstrumentmodels.jl")
end

@testset "Storage Utilities" begin
    include("teststorageutils.jl")
end

@testset "Utility Functions" begin
    include("testutilities.jl")
end
