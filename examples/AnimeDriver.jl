ENV["JULIA_CONDAPKG_BACKEND"] = "Null" # never install Conda packages; run script after activating the correct Conda environment externally

using ArgParse
using YAML
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
            help = "Overwrite existing arguments"
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

# change to output directory, optionally creating it if it does not exist
isdir(args["outdir"]) || mkdir(args["outdir"])
@info("Changing working directory to $(args["outdir"])")
cd(args["outdir"])

# create a new empty MS -- check if an MS of the same name exists and if yes, delete before creation
isdir(yamlconf["msname"]) ? (args["clobber"] || error("$(yamlconf["msname"]) exists! Not overwriting.")) : run(`rm -rf $(yamlconf["msname"])`)
generatems(yamlconf, ",", false, template) # comma-separated; do not ignore repeated delimiters

# call wscean to predict visibilities -- TODO: read in hdf5 model and create fitsdir with sky models
if yamlconf["skymodelmode"] == "fits"
    runwsclean(yamlconf["msname"], yamlconf["fits"], yamlconf["polarized"], yamlconf["channelgroups"], yamlconf["osfactor"])
elseif yamlconf["skymodelmode"] == "hdf5"
    runwsclean(yamlconf["msname"], yamlconf["hdf5"], yamlconf["polarized"], yamlconf["channelgroups"], yamlconf["osfactor"])
else
    error("Unrecognised value \"$(yamlconf["skymodelmode"])\" for modelmode! Allowed values are 'hdf5' or 'fits'.")
end

# load ms data into custom struct
#observation = loadobs(yamlconf["msname"], yamlconf["stations"])
#@info("Measurement Set and station info loaded into memory for processing")

# add thermal noise
#yamlconf["thermalnoise"]["enable"] && thermalnoise(observation, yamlconf)

# Change back to original working directory
@info("Changing working directory back to $(startdir)")
cd(startdir)
@info("Anime run completed successfully 💯")
