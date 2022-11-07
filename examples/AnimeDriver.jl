ENV["JULIA_CONDAPKG_BACKEND"] = "Null" # never install Conda packages; run script after activating the correct Conda environment externally

using ArgParse
using JSON3
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


# Change to output directory -- all operations are performed within the output directory
startdir = pwd()
isdir(args["outdir"]) || mkdir(args["outdir"])
@info("Changing working directory to $(args["outdir"])")
cd(args["outdir"])
# ensure argparse paths are absolute
config = abspath(startdir, args["config"])
template = abspath(startdir, args["template"])

# load input json file
json_string = read(config, String)
jsonpars = JSON3.read(json_string)

# create an empty MS
genms(jsonpars, template, args["clobber"])

# call wscean to predict visibilities
runwsclean(jsonpars.msname, jsonpars.fitsdir, Bool(jsonpars.toggle_polmodel), jsonpars.channelgroups, jsonpars.osfactor)

# Change back to original working directory
@info("Changing working directory back to $(startdir)")
cd(startdir)
