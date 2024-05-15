# ### A note on generating tutorial outputs
# All the tutorial outputs along with the entire documentation can be generated locally from the source code as shown below:
# ```bash
# $ git clone https://github.com/iniyannatarajan/Anime.jl.git
# $ cd Anime.jl
# $ julia --project=docs/
# (docs) pkg> instantiate
# (docs) pkg> dev .
# julia> include("docs/make.jl")
# ```
# The outputs are generated in the `docs/build` directory.

# # Creating data sets

# `Anime` uses the `CASA` Measurement Set (MS) data format as on-disk storage format from which data are read and written back to.
# Conversion to and from the more traditional uvfits format is also supported. While there is no standard modern data format yet
# for storing VLBI data, support for a handful of commonly used data formats will be made available.

# Assuming `Anime` is installed and the necessary dependencies are met, the following code snippets demonstrate how to create an MS from scratch and from an existing uvfits file.
# We load `Anime` first.
# ```julia
# using Anime
# ```

relativepath = dirname(dirname(dirname(@__DIR__))) # hide

include(joinpath(relativepath, "src", "Anime.jl")) # hide
using .Anime # hide

# ## In manual mode

# In manual mode, the function [`createmsfromconfig`](@ref Anime.createmsfromconfig) is used to create an MS from scratch,
# with the observation and site parameters passed as arguments. The `casatools` python library is used under the hood.
# In the following, replace `relativepath` with the path to the source code of `Anime.jl` which contains sample input files
# under `inputs/`.
obsconfig = Dict()
obsconfig["msname"] = "eht.ms"
obsconfig["mode"] = "manual"
obsconfig["stations"] = joinpath(relativepath, "inputs", "eht_2017.stations")
obsconfig["casaanttemplate"] = joinpath(relativepath, "inputs", "antenna_table.template")
obsconfig["spw"] = Dict()
obsconfig["spw"]["centrefreq"] = [229.0e9]
obsconfig["spw"]["bandwidth"] = [2.0e9]
obsconfig["spw"]["channels"] = [32]
obsconfig["source"] = Dict{String, Any}("M87" => Dict{String, String}("RA"=>"12h30m49.42", "Dec"=>"+12.23.28.04", "epoch"=>"J2000"))
obsconfig["starttime"] = "UTC,2021/04/28/00:00:00.00"
obsconfig["exposure"] = 1.0
obsconfig["scans"] = 2
obsconfig["scanlengths"] = [900.0, 600.0]
obsconfig["scanlag"] = 300.0
obsconfig["autocorr"] = false
obsconfig["telescopename"] = "VLBA"
obsconfig["feed"] = "perfect R L"
obsconfig["shadowlimit"] = 1e-6
obsconfig["elevationlimit"] = "10deg"
obsconfig["stokes"] = "RR RL LR LL"

delim = ","
ignorerepeated = false

stationinfo = readstationinfo(obsconfig["stations"], delim=",", ignorerepeated=false) # read station info file

createmsfromconfig(obsconfig, stationinfo)

# The above step creates a fully functional MS that can be used for further processing.

rm("eht.ms", force=true, recursive=true) # hide

rm("ANTENNA_eht_2017", force=true, recursive=true) # hide


# ## In uvfits mode

# This method accepts an existing uvfits file (e.g. output by `eht-imaging`) and uses `CASA` to convert between the two formats.
# This is done via [`createmsfromuvfits`](@ref Anime.createmsfromuvfits).
# As before, replace `relativepath` with the path to the source code of `Anime.jl` which contains sample input files
# under `inputs/`:
uvfits = joinpath(relativepath, "inputs", "uvfitsfiles", "hops_lo_3601_M87+zbl-dtcal_selfcal.uvfits")
msname = "eht.ms"
mode = "uvfits"

createmsfromuvfits(uvfits, msname, mode)

# It is the responsibility of the user to ensure that the input uvfits file contains all the necessary information that `CASA` would need to create an MS.
# `eht-imaging` output files are consistent with these specifications.

# A helper function to convert an MS back to uvfits format is also provided:
msname = "eht.ms" # the MS we just created above
uvfits = "eht.uvfits"
datacolumn = "data"

createuvfitsfromms(msname, uvfits, datacolumn)

rm("eht.ms", force=true, recursive=true) # hide

rm("eht.uvfits", force=true, recursive=true) # hide

rm("ANTENNA_eht_2017", force=true, recursive=true) # hide
