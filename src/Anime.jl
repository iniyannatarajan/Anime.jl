module Anime

using PythonCall
using Logging

include("generatems/genms.jl")
include("coherency/wscleanpredict.jl")

export genms, runwsclean

end
