#ENV["JULIA_CONDAPKG_BACKEND"] = "Null" # uncomment to never install Conda packages if an external custom Conda environment is desired

using ArgParse
using Logging
using HDF5

include("../src/Anime.jl")
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
Anime.Observe.generatems(config, delim=",", ignorerepeated=false) # comma-separated; do not ignore repeated delimiters

# call wscean to predict visibilities
Anime.Observe.computecoherency(config)

# load ms data into custom struct
obs = Anime.Observe.loadobs(config, delim=",", ignorerepeated=false)

# add corruptions
#addcorruptions(obs)
# create HDF5 file to store all corruptions
@info("Initialising empty HDF5 file to store propagation path effects...")
fid = h5open(obs.yamlconf["hdf5corruptions"], "w") # using mode "w" to destroy existing contents
close(fid)

# add tropospheric effects
obs.yamlconf["troposphere"]["enable"] && Anime.ModelInstrument.troposphere(obs)

# add instrumental polarization
obs.yamlconf["instrumentalpol"]["enable"] && Anime.ModelInstrument.instrumentalpol(obs)

# add pointing errors
obs.yamlconf["pointing"]["enable"] && Anime.ModelInstrument.pointing(obs)

# add station gains
obs.yamlconf["stationgains"]["enable"] && Anime.ModelInstrument.stationgains(obs; draw=true)

# add bandpasses
obs.yamlconf["bandpass"]["enable"] && Anime.ModelInstrument.bandpass(obs)

# add thermal noise
obs.yamlconf["thermalnoise"]["enable"] && Anime.ModelInstrument.thermalnoise(obs)

# compute weights and write everything to disk
Anime.ModelInstrument.postprocessms(obs)

# Change back to original working directory
@info("Changing working directory back to $startdir")
cd(startdir)
@info("Anime.jl observation completed successfully 📡")