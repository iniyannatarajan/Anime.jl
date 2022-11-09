module Anime

using YAML
using Logging
using PythonCall

include("observation/observation.jl")
include("corruption/corruption.jl")

end
