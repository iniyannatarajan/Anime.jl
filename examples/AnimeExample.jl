#ENV["JULIA_CONDAPKG_BACKEND"] = "Null" # uncomment to never install Conda packages if an external custom Conda environment is desired

using ArgParse
using Logging
using HDF5
using YAML

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

y = YAML.load_file(config, dicttype=Dict{String,Any})

# generate new MS
if y["mode"] == "manual"
    msfromconfig(y["msname"], y["mode"], y["stations"], y["casaanttemplate"], y["spw"]["centrefreq"], y["spw"]["bandwidth"], y["spw"]["channels"],
    y["source"], y["starttime"], y["exposure"], y["scans"], y["scanlengths"], y["scanlag"]; autocorr=y["autocorr"], telescopename=y["telescopename"],
    feed=y["feed"], shadowlimit=y["shadowlimit"], elevationlimit=y["elevationlimit"], stokes=y["stokes"], delim=",", ignorerepeated=false)

elseif y["mode"] == "uvfits"
    msfromuvfits(y["uvfits"], y["msname"], y["mode"], y["stations"], delim=",", ignorerepeated=false)
else
    error("MS generation mode '$(y["mode"])' not recognised 🤷")
end

# call wscean to compute source coherency
run_wsclean(y["msname"], y["skymodel"], y["polarized"], y["channelgroups"], y["osfactor"])

# load ms data into custom struct
obs = loadms(y, y["msname"], y["stations"], Int(y["corruptseed"]), Int(y["troposphere"]["tropseed"]), delim=",", ignorerepeated=false)

# make diagnostic plots of uncorrupted data
y["diagnostics"] && plotvis(obs, saveprefix="modelvis_") #plotvis(obs.data, obs.flag, obs.uvw, obs.chanfreqvec, obs.numchan, saveprefix="beforepropagation_")

# add corruptions
#addcorruptions(obs)

h5file = "insmodel.h5"

# add tropospheric effects
# TODO must pass h5file name to troposphere for now
y["troposphere"]["enable"] && troposphere(obs, h5file=h5file)

# add instrumental polarization
y["instrumentalpol"]["enable"] && instrumentalpol(obs, h5file=h5file)

# add pointing errors
if y["pointing"]["enable"]
    #pointing(obs, h5file=h5file)
    pointing(obs.stationinfo, obs.scanno, obs.chanfreqvec, y["pointing"]["interval"], y["pointing"]["mode"], obs.exposure, obs.times, obs.rngcorrupt,
    obs.antenna1, obs.antenna2, obs.data, obs.numchan, h5file=h5file)
    y["diagnostics"] && plotpointingerrors(obs)
end

# add station gains
if y["stationgains"]["enable"]
    stationgains(obs.scanno, obs.times, obs.exposure, obs.data, obs.stationinfo, y["stationgains"]["mode"],
    obs.rngcorrupt, obs.antenna1, obs.antenna2, obs.numchan, h5file=h5file)
    y["diagnostics"] && plotstationgains(obs)
end

# add bandpasses
if y["bandpass"]["enable"]
    bandpass(y["bandpass"]["bandpassfile"], obs.data, obs.stationinfo, obs.rngcorrupt, obs.antenna1, obs.antenna2,
    obs.numchan, obs.chanfreqvec, h5file=h5file)
    y["diagnostics"] && plotbandpass(obs)
end

# add thermal noise
y["thermalnoise"]["enable"] && thermalnoise(obs.times, obs.antenna1, obs.antenna2, obs.data,
y["correff"], obs.exposure, obs.chanwidth, obs.rngcorrupt, obs.stationinfo.sefd_Jy, h5file=h5file)

# make diagnostic plots
y["diagnostics"] && plotvis(obs, saveprefix="datavis_") #plotvis(obs.data, obs.flag, obs.uvw, obs.chanfreqvec, obs.numchan, saveprefix="afterpropagation_")

# compute weights and write everything to disk
postprocessms(obs, h5file=h5file)

# Change back to original working directory
@info("Changing working directory back to $startdir")
cd(startdir)
@info("Anime.jl observation completed successfully 📡")