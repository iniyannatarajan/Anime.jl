module Observe

using CSV
using YAML
using Random
using DataFrames
using Casacore.Tables: Tables as CCTables, Table as CCTable
using PythonCall

importuvfits = pyimport("casatasks" => "importuvfits")
table = pyimport("casatools" => "table")
measures = pyimport("casatools" => "measures")
simulator = pyimport("casatools" => "simulator")
tb = table()
me = measures()
sm = simulator()

include(joinpath(@__DIR__, "generatems.jl"))
include(joinpath(@__DIR__, "coherency.jl"))
include(joinpath(@__DIR__, "loadobs.jl"))

end