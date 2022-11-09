export thermalnoise, skynoise

using CSV
using DataFrames
using Casacore.Tables: Tables, Table

function thermalnoise(stations::String, delim::String, ignorerepeated::Bool, correff::Float64)
    """
    Add thermal noise to visibilities
    """
    df = CSV.read(stations, DataFrame; delim=delim, ignorerepeated=ignorerepeated)

end

function skynoise()
end
