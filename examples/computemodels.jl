# # Compute instrument models

# The primary goal of `Anime` is to generate instrument models tailored to observation specifications. The models are optionally stored in HDF5 format
# and plotted. Plotting functions using `Plots.jl` are provided.

# All functions that follow require loading the following modules:
relativepath = "../../../"

include(joinpath(relativepath, "src", "Anime.jl"))
using .Anime
using HDF5
using Plots

# We also load an *observation* from an existing MS to operate on.
msname = joinpath(relativepath, "test", "data", "ehtuvf.ms")
stations = joinpath(relativepath, "inputs", "eht_2017.stations")
corruptseed = 456
tropseed = 54364
tropwetonly =false
correff = 0.88
tropattenuate = true
tropskynoise = true
tropmeandelays = true
tropturbulence = true
polvisframe = "sky"
polmode = "gp"
ptginterval = 5.0
ptgscale = 2.0
ptgmode = "gp"
gainsmode = "gp"
bpfile = joinpath(relativepath, "test", "data", "eht_2017.bandpass")
delim = ","
ignorerepeated = false

obs = loadms(msname, stations, corruptseed, tropseed, tropwetonly, correff, tropattenuate, tropskynoise, tropmeandelays, tropturbulence, polvisframe,
polmode, ptginterval, ptgscale, ptgmode, gainsmode, bpfile, delim=",", ignorerepeated=false)

# ## Atmospheric models
# At mm-wavelengths (230 GHz), the troposphere has significant effects on signal propagation. `Anime` re-implements in `Julia` all the tropospheric effects
# simulated by `MEQSv2`[^1]. The advantage here is that apart from being faster and not requiring a regular grid of complex visibilities, they can be
# called from other imaging and calibration packages when necessary.

# The function [`troposphere!`](@ref Anime.troposphere!) computes various tropospheric effects based on the flags set when `loadms` is called.

# !!! note
#     The external program `AATM` is requred to generate quantities related to tropospheric absorption and dispersion.
#     Here we pass pre-existing CSV output from `AATM` to `absorptionfile` and `dispersivefile`. The elevation angles of all antennas
#     during the course of the observation can also be preloaded when `casatools` are unavailable to generate them from scratch.

h5file = "sample.h5"
absorptionfile = joinpath(relativepath, "test", "data", "absorption1.csv")
dispersivefile = joinpath(relativepath, "test", "data", "dispersive1.csv")
elevfile = joinpath(relativepath, "test", "data", "insmodeluvf.h5")

troposphere!(obs, h5file, absorptionfile=absorptionfile, dispersivefile=dispersivefile, elevfile=elevfile)

# This computes the delays introduced by the mean and turbulent components of the troposphere, along with attenuation due to opacity and increase in
# system noise.

# `Anime` provides plotting functions that help visualize these models. In the following, the gaps in the plotted curves signify lags between two
# observing scans.

# For example, to plot the elevation angles by station we can just do
plotelevationangle(elevfile, obs.scanno, obs.times, obs.stationinfo.station, save=false)

# The transmission values computed can be plotted using
plottransmission(h5file, obs.stationinfo.station, obs.times, obs.chanfreqvec, save=false)
# Since this is a channel-averaged data set, the frequency-dependent transmission reduces to a single curve per station.

# The delays due to the mean component of the troposphere can be plotted as follows:
plotmeandelays(h5file, obs.stationinfo.station, obs.times, obs.chanfreqvec, save=false)
#-
rm(h5file) # hide
rm("atm.csv") # hide

# ## Primary beam response
# We model the effects of antenna pointing offsets caused by various mechanical or electronic effects. The size of the primary beam affects how much the errors
# in antenna pointing attenuate the complex visibilities. `Anime` models the pointing offsets as Gaussian processes which are useful for modelling smooth
# variations in pointing over time.

# The function [`pointing!`](@ref Anime.pointing!) is used to compute and apply pointing models to data.
pointing!(obs, h5file=h5file)
# Note that this method is a shorthand for another method with multiple arguments that provides more fine-grained control over the input parameters.

# We now plot the pointing model generated:
plotpointingerrors(h5file, obs.scanno, obs.stationinfo.station, save=false)
# Mispointings of station LM (Large Millimeter Telescope, Mexico), the largest dish in the array, result in the largest attenuation of amplitude.
rm(h5file) # hide

# ### References
# [^1]: Natarajan I. et al. MeqSilhouette v2 (2022) [MNRAS](https://academic.oup.com/mnras/article/512/1/490/6537429)