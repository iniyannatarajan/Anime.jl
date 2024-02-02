# # Read/write MS

# `Anime` stores all information necessary to compute and apply instrument models to data in a custom structure of Arrays
# (see [`CjlObservation`](@ref Anime.CjlObservation)). All relevant data are read from the MS into this struct using [`loadms`](@ref Anime.loadms).

# Load necessary modules

relativepath = "../../../"

include(joinpath(relativepath, "src", "Anime.jl"))
using .Anime

# For this example, we will load an MS from the `test/data/` directory. 
msname = joinpath(relativepath, "test", "data", "eht1.ms")
stations = joinpath(relativepath, "inputs", "eht_2017.stations")
corruptseed = 456
tropseed = 54364
tropwetonly = false
correff = 0.88
tropattenuate = true
tropskynoise = true
tropmeandelays = true
tropturbulence = true
polvisframe = "sky"
polmode = "gp"
ptginterval = 5.0
ptgscale = 2.0
ptgmode = "gp"
gainsmode = "gp"
bpfile = joinpath(relativepath, "test", "data", "eht_2017.bandpass")
delim = ","
ignorerepeated = false

obs = loadms(msname, stations, corruptseed, tropseed, tropwetonly, correff, tropattenuate, tropskynoise, tropmeandelays, tropturbulence, polvisframe,
polmode, ptginterval, ptgscale, ptgmode, gainsmode, bpfile, delim=",", ignorerepeated=false)

# All subsequent computations are performed on this structure until the user calls [`postprocessms`](@ref Anime.postprocessms) 
# to write the results back to the MS. This function optionally accepts an HDF5 file containing the noise to be used for populating
# the WEIGHT and SIGMA arrays in the MS. If no HDF5 file is provided, these arrays are initialized to 1.0 and 0.0 respectively.

# We first load a pre-existing HDF5 file containing instrument models generated according to the dimensions required by `eht1.ms` that was loaded 
# in the previous step.

h5file = joinpath(relativepath, "test", "data", "insmodel1.h5")

postprocessms(obs, h5file=h5file)

# The above function call writes the contents of `obs` back into the MS.

# !!! note
#     Currently `Casacore.jl` does not provide a way to cleanly close and release the read-write lock to an MS. If further processing is required using a different library such as `casatools` it is recommended to do that in a new REPL session.