using Anime
using Test

#=@testset "Create MS" begin
    include("testcreatems.jl")
end=#

#=@testset "Compute Coherency" begin
    include("testcomputecoherency.jl")
end=#

@testset "Load Observation" begin
    include("testloadobservation.jl")
end

#=@testset "Model Instrument" begin
    #include("testthermalnoise.jl")
    #include("teststationgains.jl")
    #include("testbandpass.jl")
    #include("testbeam.jl")
    #include("testinstrumentalpol.jl")
    #include("testtroposphere.jl")
end=#