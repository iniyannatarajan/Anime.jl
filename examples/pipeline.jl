# # Running in pipeline mode

# `Anime` can be run in pipeline mode with little to no user interaction once the pipeline has started execution. This mode can include everything from
# creating a new data set from scratch to computing coherency matrices, computing and applying instrumental models to complex visibilities, and writing
# gain tables and making diagnostic plots. Note that some steps may involve the use of external software such as `WSClean`.

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

# Now we load the YAML file containing the observation and instrument modelling parameters. While this is not mandatory, it is far easier to keep track of the
# input settings when a configuration file is used.
# ```julia
# y = YAML.load_file(config, dicttype=Dict{String,Any}) # load YAML config file
# h5file = "insmodel.h5"
# ```

# Generate a new data set from either the YAML file observation parameters or a previously existing UVFITS file. An ASCII file containing station information
# ([Site-and-station-parameters](@ref)) is necessary in both cases. Also see [Creating-data-sets](@ref).
# ```julia
# if y["mode"] == "manual"
#     msfromconfig(y["msname"], y["mode"], y["stations"], y["casaanttemplate"], y["spw"]["centrefreq"], y["spw"]["bandwidth"], y["spw"]["channels"],
#     y["source"], y["starttime"], y["exposure"], y["scans"], y["scanlengths"], y["scanlag"]; autocorr=y["autocorr"], telescopename=y["telescopename"],
#     feed=y["feed"], shadowlimit=y["shadowlimit"], elevationlimit=y["elevationlimit"], stokes=y["stokes"], delim=",", ignorerepeated=false)
# elseif y["mode"] == "uvfits"
#     msfromuvfits(y["uvfits"], y["msname"], y["mode"])
# else
#     error("MS generation mode '$(y["mode"])' not recognised 🤷")
# end
# ``` 

# Once a new MS has been generated, we can populate the source coherency in the `DATA` column using any visibility prediction software. Here, we use `WSClean`
# to compute source coherency (also see [Compute-coherency-matrix](@ref)):
# ```julia
# run_wsclean(y["msname"], y["skymodel"], y["polarized"], y["channelgroups"], y["osfactor"])
# ```

# Now we load the data into a custom struct:
# ```julia
# obs = loadms(y["msname"], y["stations"], Int(y["corruptseed"]), Int(y["troposphere"]["tropseed"]), y["troposphere"]["wetonly"], y["correff"],
# y["troposphere"]["attenuate"], y["troposphere"]["skynoise"], y["troposphere"]["meandelays"], y["troposphere"]["turbulence"],
# y["instrumentalpolarization"]["visibilityframe"], y["instrumentalpolarization"]["mode"], y["pointing"]["interval"], y["pointing"]["scale"], 
# y["pointing"]["mode"], y["stationgains"]["mode"], y["bandpass"]["bandpassfile"], delim=",", ignorerepeated=false)
# ```
 
# We can now start adding the individual instrument models. For more details, see [Compute-instrument-models](@ref).
# ```julia
# y["troposphere"]["enable"] && troposphere!(obs, h5file)
# y["instrumentalpolarization"]["enable"] && instrumentalpolarization!(obs, h5file=h5file)
# y["pointing"]["enable"] && pointing!(obs, h5file=h5file)
# y["stationgains"]["enable"] && stationgains!(obs, h5file=h5file)
# y["bandpass"]["enable"] && bandpass!(obs, h5file=h5file)
# y["thermalnoise"]["enable"] && thermalnoise!(obs, h5file=h5file)
# ```

# Compute weights, write everything back to disk, and convert to uvfits
# ```julia
# postprocessms(obs, h5file=h5file)
# mstouvfits(y["msname"], "test.uvfits", "corrected")
# ```

# Change back to original working directory and exit
# ```julia
# cd(startdir)
# @info("Anime.jl observation completed successfully 📡")
# ```
