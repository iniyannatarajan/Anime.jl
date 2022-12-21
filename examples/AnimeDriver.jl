ENV["JULIA_CONDAPKG_BACKEND"] = "Null" # never install Conda packages; run script after activating the correct Conda environment externally

using ArgParse
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
        "outdir"
            help = "Directory to hold output products"
            required = true
        "--clobber", "-c"
            action = :store_true
            help = "Delete and create output directory anew"
    end

    return parse_args(s)
end

# create parser
args = create_parser()

# change working directory to the user-specified output directory
startdir = pwd() # store original working directory
config = abspath(startdir, args["config"])
outdir = abspath(startdir, args["outdir"])

# create a new empty output directory
if isdir(outdir)
    if args["clobber"]
        run(`rm -rf $(outdir)`)
	mkdir(outdir)
    else 
	error("$outdir exists but -c option is not given 🤷") 
    end
else
    mkdir(outdir)
end
@info("Changing working directory to $outdir")
cd(outdir)

# create new ms
@time generatems(config, delim=",", ignorerepeated=false) # comma-separated; do not ignore repeated delimiters

# call wscean to predict visibilities
@time predict_visibilities(config)

# load ms data into custom struct
@time observation = loadobs(config, delim=",", ignorerepeated=false)

# add corruptions
addcorruptions(observation)

# Change back to original working directory
@info("Changing working directory back to $startdir")
cd(startdir)
@info("Anime.jl observation completed successfully 📡")