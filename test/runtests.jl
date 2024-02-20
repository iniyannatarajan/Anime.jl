include("../src/Anime.jl")
using .Anime
#using Anime
using Test

using CairoMakie, ColorSchemes, LaTeXStrings
using Tables
using Logging
using CSV
using Glob
using HDF5
using YAML
using Random
using Statistics
using SpecialFunctions
using DataFrames
using Distributions
using DocStringExtensions
using LinearAlgebra
using Casacore.Tables: Tables as CCTables, Table as CCTable
using PythonCall

quanta = pyimport("casatools" => "quanta")
importuvfits = pyimport("casatasks" => "importuvfits")
exportuvfits = pyimport("casatasks" => "exportuvfits")
table = pyimport("casatools" => "table")
measures = pyimport("casatools" => "measures")
simulator = pyimport("casatools" => "simulator")
qa = quanta()
tb = table()
me = measures()
sm = simulator()

@testset "Create data set" begin
    include("testcreatems.jl")
end

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