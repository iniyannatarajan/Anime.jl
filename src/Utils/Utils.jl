module Utils

using Random
using Casacore.Tables: Tables as CCTables, Table as CCTable
using PythonCall

quanta = pyimport("casatools" => "quanta")
measures = pyimport("casatools" => "measures")
qa = quanta()
me = measures()

using ..Observe: CjlObservation

include(joinpath(@__DIR__, "util.jl"))

end