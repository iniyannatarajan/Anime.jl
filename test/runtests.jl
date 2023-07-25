using Anime
using Test
using YAML
using Random

#=@testset "Create MS" begin
    include("testcreatems.jl")
end=#

#=@testset "Compute Coherency" begin
    include("testcomputecoherency.jl")
end=#

@testset "Load Observation" begin
    include("testloadobs.jl")
end

@testset "Utility Functions" begin
    include("testutilities.jl")
end
