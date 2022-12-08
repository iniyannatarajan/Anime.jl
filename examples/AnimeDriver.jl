ENV["JULIA_CONDAPKG_BACKEND"] = "Null" # never install Conda packages; run script after activating the correct Conda environment externally

using ArgParse
using YAML
using Random
using Logging

include("../src/Anime.jl")
using .Anime

# create argument parser
function create_parser()
    s = ArgParseSettings()

    @add_arg_table s begin
        "config"
            help = "Input JSON file name with observation settings"
            required = true
        "--template", "-t"
            help = "Input CASA ANTENNA table to use as template"
            arg_type = String
        "--outdir", "-o"
            help = "Directory to hold output products"
            arg_type = String
        "--clobber", "-c"
            action = :store_true
            help = "Delete and create output directory anew"
    end

    return parse_args(s)
end

# create parser
args = create_parser()

# get absolute paths of the input files
startdir = pwd()
config = abspath(startdir, args["config"])
template = abspath(startdir, args["template"])

# load input YAML file
yamlconf = YAML.load_file(config, dicttype=Dict{String,Any})

# create a new empty output dir if clobber==true
if isdir(args["outdir"])
    if args["clobber"]
        run(`rm -rf $(args["outdir"])`)
	mkdir(args["outdir"])
    else 
	error("Output dir '$(args["outdir"])' exists but clobber=false 🤷") 
    end
else
    mkdir(args["outdir"])
end

@info("Changing working directory to $(args["outdir"])")
cd(args["outdir"])

# create new ms
@time generatems(yamlconf, ",", false, template) # comma-separated; do not ignore repeated delimiters

# call wscean to predict visibilities -- TODO: read in hdf5 model and create fitsdir with sky models
if yamlconf["skymodelmode"] == "fits"
    @time runwsclean(yamlconf["msname"], yamlconf["fitssky"], yamlconf["polarized"], yamlconf["channelgroups"], yamlconf["osfactor"])
elseif yamlconf["skymodelmode"] == "hdf5"
    @time runwsclean(yamlconf["msname"], yamlconf["hdf5sky"], yamlconf["polarized"], yamlconf["channelgroups"], yamlconf["osfactor"])
else
    error("Unrecognised value \"$(yamlconf["skymodelmode"])\" for modelmode! Allowed values are 'hdf5' or 'fits'.")
end

# load ms data into custom struct
#observation, stationinfo = loadobs(yamlconf["msname"], yamlconf["stations"], ",", false)
@time observation = loadobs(yamlconf, ",", false)

# add corruptions
addcorruptions(observation)

# Change back to original working directory
@info("Changing working directory back to $(startdir)")
cd(startdir)
@info("Anime observation completed successfully 📡")
