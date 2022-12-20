ENV["JULIA_CONDAPKG_BACKEND"] = "Null" # never install Conda packages; run script after activating the correct Conda environment externally

using ArgParse
using YAML
using Random
using Logging

@info("Including Anime.jl ...")
@time include("../src/Anime.jl")
using .Anime

# create argument parser
function create_parser()
    s = ArgParseSettings()

    @add_arg_table s begin
        "config"
            help = "Input YAML file name with observation configuration"
            required = true
        "--clobber", "-c"
            action = :store_true
            help = "Delete and create output directory anew"
    end

    return parse_args(s)
end

# create parser
args = create_parser()

startdir = pwd() # store original working directory
config = abspath(startdir, args["config"])
yamlconf = YAML.load_file(config, dicttype=Dict{String,Any}) # load input YAML file

template = abspath(startdir, yamlconf["casatemplate"]) # get absolute path of casa antenna template

# create a new empty output dir if clobber==true
if isdir(yamlconf["outdir"])
    if args["clobber"]
        run(`rm -rf $(yamlconf["outdir"])`)
	mkdir(yamlconf["outdir"])
    else 
	error("Output dir '$(yamlconf["outdir"])' exists but clobber=false 🤷") 
    end
else
    mkdir(yamlconf["outdir"])
end

@info("Changing working directory to $(yamlconf["outdir"])")
cd(yamlconf["outdir"])

# create new ms
@time generatems(yamlconf, ",", false, template) # comma-separated; do not ignore repeated delimiters

# call wscean to predict visibilities
@time predict_visibilities(yamlconf)

# load ms data into custom struct
#observation, stationinfo = loadobs(yamlconf["msname"], yamlconf["stations"], ",", false)
@time observation = loadobs(yamlconf, delim=",", ignorerepeated=false)

# add corruptions -- TODO currently, the stationinfo file and the MS ANTENNA table should be in the same order! Make this more flexible!!!
addcorruptions(observation)

# Change back to original working directory
@info("Changing working directory back to $(startdir)")
cd(startdir)
@info("Anime.jl observation completed successfully 📡")