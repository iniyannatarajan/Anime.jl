# # Creating data sets

# `Anime` uses the `CASA` Measurement Set (MS) data format as on-disk storage format from which data are read and written back to.
# Conversion to and from the more traditional uvfits format is also supported. While there is no standard format for storing VLBI data, support for a handful of
# commonly used data formats will be made available.

# We load the necessary modules first.

include("../../../src/Anime.jl")
using .Anime
#using Anime
using YAML

# Now load the input config file that contains the observation parameters. Sample files are provided in the source code under the `inputs/` directory.
config = "../../../inputs/config.yaml"
y = YAML.load_file(config, dicttype=Dict{String,Any})

# ## In manual mode

# The first method obtains all the observation and site parameters from a `YAML` file and uses `CASA` to build
# an MS. To generate a new MS from scratch we call [`msfromconfig`](@ref Anime.msfromconfig):
msfromconfig("eht.ms", "manual", y["stations"], y["casaanttemplate"], y["spw"]["centrefreq"], y["spw"]["bandwidth"], y["spw"]["channels"],
y["source"], y["starttime"], y["exposure"], y["scans"], y["scanlengths"], y["scanlag"]; autocorr=y["autocorr"], telescopename=y["telescopename"],
feed=y["feed"], shadowlimit=y["shadowlimit"], elevationlimit=y["elevationlimit"], stokes=y["stokes"], delim=",", ignorerepeated=false)

# The above step creates a fully functional MS that can be used for further processing.

rm("eht.ms", force=true, recursive=true) # hide

# ## In uvfits mode

# The second method accepts an existing uvfits file (e.g. output by `eht-imaging`) and uses `CASA` to convert between the two formats.
# This is done via [`msfromuvfits`](@ref Anime.msfromuvfits):

msfromuvfits(y["uvfits"], "eht.ms", "uvfits", y["stations"], delim=",", ignorerepeated=false)

# It is the responsibility of the user to ensure that the input uvfits file contains all the necessary information that `CASA` would expect.
# `eht-imaging` output files are consistent with these specifications.

rm("eht.ms", force=true, recursive=true) # hide

