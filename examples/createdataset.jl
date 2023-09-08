# # Creating data sets

# `Anime` uses the `CASA` Measurement Set (MS) data format as on-disk storage format from which data are read and written back to.
# Conversion to and from the more traditional uvfits format is also supported. While there is no standard format for storing VLBI data, support for a handful of
# commonly used data formats will be made available.

# We load the necessary modules first.

relativepath = "../../../"

include(joinpath(relativepath, "src", "Anime.jl"))
using .Anime
#using Anime

# ## In manual mode

# In manual mode, the function [`msfromconfig`](@ref Anime.msfromconfig) is used to create an MS from scratch, with the observation and site parameters
# passed as arguments. The `casatools` python library is used under the hood.
msname = "eht.ms"
mode = "manual"
stations = joinpath(relativepath, "inputs", "eht_2017.stations")
casaanttemplate = joinpath(relativepath, "inputs", "antenna_table.template")
spw_centrefreq = [229.0e9]
spw_bw = [2.0e9]
spw_channels = [32]
sourcedict = Dict{String, Any}("M87" => Dict{String, String}("RA"=>"12h30m49.42", "Dec"=>"+12.23.28.04", "epoch"=>"J2000"))
starttime = "UTC,2021/04/28/00:00:00.00"
exposure = 1.0
scans = 2
scanlengths = [900.0, 600.0]
scanlag = 300.0
autocorr = false
telescopename = "VLBA"
feed = "perfect R L"
shadowlimit = 1e-6
elevationlimit = "10deg"
stokes = "RR RL LR LL"
delim = ","
ignorerepeated = false

msfromconfig(msname, mode, stations, casaanttemplate, spw_centrefreq, spw_bw, spw_channels, sourcedict, starttime, exposure, scans, scanlengths, scanlag;
 autocorr=autocorr, telescopename=telescopename, feed=feed, shadowlimit=shadowlimit, elevationlimit=elevationlimit, stokes=stokes, 
 delim=delim, ignorerepeated=ignorerepeated)

# The above step creates a fully functional MS that can be used for further processing.

rm("eht.ms", force=true, recursive=true) # hide

rm("ANTENNA_eht_2017", force=true, recursive=true) # hide


# ## In uvfits mode

# This method accepts an existing uvfits file (e.g. output by `eht-imaging`) and uses `CASA` to convert between the two formats.
# This is done via [`msfromuvfits`](@ref Anime.msfromuvfits):
uvfits = joinpath(relativepath, "inputs", "uvfitsfiles", "hops_lo_3601_M87+zbl-dtcal_selfcal.uvfits")
msname = "eht.ms"
mode = "uvfits"

msfromuvfits(uvfits, msname, mode)

# It is the responsibility of the user to ensure that the input uvfits file contains all the necessary information that `CASA` would need to create an MS.
# `eht-imaging` output files are consistent with these specifications.

# A helper function to convert an MS back to uvfits format is also provided:
msname = "eht.ms"
uvfits = "eht.uvfits"
datacolumn = "data"

mstouvfits(msname, uvfits, datacolumn)

rm("eht.ms", force=true, recursive=true) # hide

rm("eht.uvfits", force=true, recursive=true) # hide

rm("ANTENNA_eht_2017", force=true, recursive=true) # hide
