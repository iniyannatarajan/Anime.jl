# SPDX-License-Identifier: MIT

# Anime - Atmospheric aNd Instrumental models in the Measurement Equation
"""Anime: Atmospheric aNd Instrumental Models in the Measurement Equation
"""
module Anime

using Tables
using Logging

include("Observe/Observe.jl")
include("Utils/Utils.jl")
include("ModelInstrument/ModelInstrument.jl")
end
