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

@testset "Stats Methods" begin
    include("teststatsutils.jl")
end

@testset "Plotting Methods" begin
    include("testplots.jl")
end