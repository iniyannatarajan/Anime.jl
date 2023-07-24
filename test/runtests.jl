using Anime
using Test

#=@testset "Create MS" begin
    include("testcreatems.jl")
end=#

#=@testset "Compute Coherency" begin
    include("testcomputecoherency.jl")
end=#

@testset "Load Observation" begin
    include("testloadobs.jl")
end

#=@testset "Utility Functions" begin
    include("testutils.jl")
end=#
