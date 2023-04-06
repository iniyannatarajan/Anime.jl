# SPDX-License-Identifier: MIT

# Anime - Atmospheric aNd Instrumental models in the Measurement Equation
"""Anime: Atmospheric aNd Instrumental Models in the Measurement Equation
"""
module Anime

using CSV
using HDF5
using YAML
using Random
using Tables
using Logging
using DataFrames
using Distributions
using Casacore.Tables: Tables as CCTables, Table as CCTable

using PythonCall
table = pyimport("casatools" => "table")
measures = pyimport("casatools" => "measures")

include("observation/observation.jl")
include("corruption/corruption.jl")

end
