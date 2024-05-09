# # Read/write MS

# `Anime` stores various arrays read from a Measurement Set to a custom struct (see [`MeasurementSet`](@ref Anime.MeasurementSet)).
# All relevant data are read from the MS into this struct using [`readms`](@ref Anime.readms).

# Load necessary modules:
# ```julia
# using Anime
# ```

relativepath = "../../../" # hide

include(joinpath(relativepath, "src", "Anime.jl")) # hide
using .Anime # hide

# For this example, we will read an MS from the `test/data/` directory. 
msname = joinpath(relativepath, "test", "data", "eht1.ms")

ms = readms(msname)

# All subsequent computations are performed on this structure until the user calls [`writems`](@ref Anime.writems) 
# to write the results back to the MS. This function optionally accepts an HDF5 file containing the noise to be used
# for populating the WEIGHT and SIGMA arrays in the MS. If no HDF5 file is provided, these arrays are initialized to
# 1.0 and 0.0 respectively.

# We first read an existing HDF5 file containing instrument models generated according to the dimensions required by `eht1.ms`
# that was read in the previous step.

h5file = joinpath(relativepath, "test", "data", "insmodel1.h5")

writems(ms, h5file=h5file)

# The above function call writes the contents of `obs` back into the MS.

# !!! note
#     Currently `Casacore.jl` does not provide a way to cleanly close and release the read-write lock to an MS. If further processing is required using a different library such as `casatools`, it is recommended to do that in a new REPL session.