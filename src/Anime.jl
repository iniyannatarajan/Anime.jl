module Anime

using CSV
using YAML
using Logging
using DataFrames
using Casacore.Tables: Tables as CasacoreTables, Table as CasacoreTable

using PythonCall
table = pyimport("casatools" => "table")
measures = pyimport("casatools" => "measures")

include("observation/observation.jl")
include("corruption/corruption.jl")

end
