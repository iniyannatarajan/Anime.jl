include("../src/Anime.jl")
using .Anime
#using Anime
#=using Test

using Plots, ColorSchemes, LaTeXStrings
using Tables
using Logging
using CSV
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
sm = simulator()=#

using BenchmarkTools

### Test instrument modelling
# load observation into custom struct
msname = joinpath(@__DIR__, "../test/data/eht1.ms")
stations = joinpath(@__DIR__, "../test/data/eht_2017.stations")
corruptseed = 456
tropseed = 873256
tropwetonly = false
correff = 0.88
tropattenuate = true
tropskynoise = true
tropmeandelays = true
tropturbulence = true
polframe = "sky"
polmode = "gp"
ptginterval = 10.0
ptgscale = 100.0
ptgmode = "gp"
stationgainsmode = "gp"
bandpassfile = joinpath(@__DIR__, "../test/data/eht_2017.bandpass")
delim = ","
ignorerepeated = false

#= @benchmark loadms($(msname), $(stations), $(corruptseed), $(tropseed), $(tropwetonly), $(correff), $(tropattenuate), $(tropskynoise), $(tropmeandelays),
$(tropturbulence), $(polframe), $(polmode), $(ptginterval), $(ptgscale), $(ptgmode), $(stationgainsmode), $(bandpassfile); delim=$(delim),
ignorerepeated=$(ignorerepeated))=#

obs = loadms(msname, stations, corruptseed, tropseed, tropwetonly, correff, tropattenuate, tropskynoise, tropmeandelays, tropturbulence, polframe, polmode,
ptginterval, ptgscale, ptgmode, stationgainsmode, bandpassfile; delim=delim, ignorerepeated=ignorerepeated)

h5file = "tropos.h5"
absorptionfile = joinpath(@__DIR__, "../test/data/absorption1.csv")
dispersivefile = joinpath(@__DIR__, "../test/data/dispersive1.csv")
elevfile = joinpath(@__DIR__, "../test/data/insmodel1.h5")

@benchmark troposphere!($(obs), $(h5file); absorptionfile=$(absorptionfile), dispersivefile=$(dispersivefile), elevfile=$(elevfile))

rm("atm.csv")
rm(h5file)