# # Generating synthetic data sets

# `Anime` can be run as a pipeline with little to no user interaction while generating synthetic data. This mode can include everything from
# creating a new data set from scratch to computing coherency matrices, computing and applying instrumental models to complex visibilities,
# and writing gain tables and making diagnostic plots.

# !!! note 
#     Some steps listed here involve the use of external software such as `WSClean` and hence this script is provided as a "passive" example that is not 
#     executed in the deployed documentation. It can be run by the user on their machine to generate synthetic data.

# Import the necessary modules
# ```julia
# using ArgParse
# using Logging
# using HDF5
# using YAML
# using Anime
# ```

# Set some convenient command-line options for the script.
# ```julia
# # create argument parser
# function create_parser()
#     s = ArgParseSettings()
# 
#     @add_arg_table s begin
#         "config"# 
#             help = "Input YAML file name with observation configuration"
#             required = true
#         "outdir"
#             help = "Directory to hold output products"
#             required = true
#         "--clobber", "-c"
#             action = :store_true
#             help = "Delete and create output directory anew"
#     end

#     return parse_args(s)
# end

# # create parser
# args = create_parser()
# ```

# Create a new working directory within which all output products will be stored.
# ```julia
# # change working directory to the user-specified output directory
# startdir = pwd() # store original working directory
# config = abspath(startdir, args["config"])
# outdir = abspath(startdir, args["outdir"])

# # create a new empty output directory
# if isdir(outdir)
#     if args["clobber"]
#         run(`rm -rf $(outdir)`)
# 	mkdir(outdir)
#     else 
# 	error("$outdir exists but -c option is not given 🤷") 
#     end
# else
#     mkdir(outdir)
# end
# @info("Changing working directory to $outdir")
# cd(outdir)
# ```

# Now we load the YAML file containing the observation and instrument modelling parameters.
# While this is not mandatory, it is far easier to keep track of the input settings when a configuration file is used.
# We also load the station information file which contains the site and station parameters.
# ```julia
# obsconfig = readobsconfig(config) # load YAML config file
# stationinfo = readstationinfo(obsconfig["stations"]) # load station info file
# h5file = "insmodel.h5"
# ```

# Generate a new data set from either the YAML file observation parameters or a previously existing UVFITS file. An ASCII file containing station information
# ([Site-and-station-parameters](@ref)) is necessary in both cases. Also see [Creating-data-sets](@ref).
# ```julia
# if obsconfig["mode"] == "manual"
#     createmsfromconfig(obsconfig["msname"], obsconfig["mode"], obsconfig["casaanttemplate"], stationinfo, obsconfig["spw"]["centrefreq"], obsconfig["spw"]["bandwidth"], obsconfig["spw"]["channels"],
#     obsconfig["source"], obsconfig["starttime"], obsconfig["exposure"], obsconfig["scans"], obsconfig["scanlengths"], obsconfig["scanlag"]; autocorr=obsconfig["autocorr"], telescopename=obsconfig["telescopename"],
#     feed=obsconfig["feed"], shadowlimit=obsconfig["shadowlimit"], elevationlimit=obsconfig["elevationlimit"], stokes=obsconfig["stokes"], delim=",", ignorerepeated=false)
# elseif obsconfig["mode"] == "uvfits"
#     createmsfromuvfits(obsconfig["uvfits"], obsconfig["msname"], obsconfig["mode"])
# else
#     error("MS generation mode '$(obsconfig["mode"])' not recognised 🤷")
# end
# ``` 

# Once a new MS has been generated, we can populate the source coherency in the `DATA` column using any visibility prediction software. Here, we use `WSClean`
# to compute source coherency (also see [Compute-coherency-matrix](@ref)):
# ```julia
# run_wsclean(obsconfig["msname"], obsconfig["skymodel"], obsconfig["polarized"], obsconfig["channelgroups"], obsconfig["osfactor"])
# ```

# Now we load the data into a custom struct:
# ```julia
# ms = readms(obsconfig["msname"])
# ```
 
# We can now start adding the individual instrument models. For more details, see [Compute-instrument-models](@ref).
# ```julia
# obsconfig["troposphere"]["enable"] && troposphere!(ms, stationinfo, obsconfig, h5file)
# obsconfig["instrumentalpolarization"]["enable"] && instrumentalpolarization!(ms, stationinfo, obsconfig, h5file=h5file)
# obsconfig["pointing"]["enable"] && pointing!(ms, stationinfo, obsconfig, h5file=h5file)
# obsconfig["stationgains"]["enable"] && stationgains!(ms, stationinfo, obsconfig, h5file=h5file)
# obsconfig["bandpass"]["enable"] && bandpassinfo = readbandpassinfo(obsconfig["bandpass"]["bandpassfile"])
# obsconfig["thermalnoise"]["enable"] && thermalnoise!(ms, stationinfo, obsconfig, h5file=h5file)
# ```

# Compute weights, write everything back to disk, and convert to uvfits. Note that the conversion to uvfits must be done after exiting the
# main pipeline script, since `Casacore` does not allow other modules (such as python `casatasks`) to access the MS in the same session.
# ```julia
# writems(ms, h5file=h5file)
# mstouvfits(obsconfig["msname"], "test.uvfits", "corrected")
# ```

# Change back to original working directory and exit
# ```julia
# cd(startdir)
# @info("Synthetic data generation completed successfully 📡")
# ```
