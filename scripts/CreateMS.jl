#ENV["JULIA_CONDAPKG_BACKEND"] = "Null" # uncomment to never install Conda packages if an external custom Conda environment is desired

using Logging
using YAML

include("../src/Anime.jl")
using .Anime

config = "../inputs/config.yaml"

y = YAML.load_file(config, dicttype=Dict{String,Any})

# generate new MS
if y["mode"] == "manual"
    msfromconfig2(y["msname"], y["mode"], y["stations"], y["casaanttemplate"], y["spw"]["centrefreq"], y["spw"]["bandwidth"], y["spw"]["channels"],
    y["source"], y["starttime"], y["exposure"], y["scans"], y["scanlengths"], y["scanlags"]; autocorr=y["autocorr"], telescopename=y["telescopename"],
    feed=y["feed"], shadowlimit=y["shadowlimit"], elevationlimit=y["elevationlimit"], stokes=y["stokes"], delim=",", ignorerepeated=false)

elseif y["mode"] == "uvfits"
    msfromuvfits(y, delim=",", ignorerepeated=false)
else
    error("MS generation mode '$(y["mode"])' not recognised 🤷")
end

@info("MS generated successfully 📡")