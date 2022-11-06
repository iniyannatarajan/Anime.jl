ENV["JULIA_CONDAPKG_BACKEND"] = "Null" # never install Conda packages; run script after activating the correct Conda environment externally

using ArgParse
using JSON3

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
        "--clobber", "-c"
            action = :store_true
            help = "Overwrite existing arguments"
    end

    return parse_args(s)
end

# create parser
args = create_parser()

# load input json file
json_string = read(args["config"], String)
jsonpars = JSON3.read(json_string)

# create an empty MS
genms(jsonpars, args["template"], args["clobber"])

# call wscean to predict visibilities
#slicewscleanoutputs(jsonpars.fitsdir, jsonpars.toggle_polmodel, jsonpars.channelgroups)
