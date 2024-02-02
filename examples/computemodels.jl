# # Compute instrument models

# The primary goal of `Anime` is to generate instrument models tailored to specific observations. The models are optionally stored in HDF5 format.
# Basic plotting functions are also provided to visualize the models. The following example demonstrates how to compute and visualize the models.

# All functions that follow require loading the following modules:
relativepath = "../../../"

include(joinpath(relativepath, "src", "Anime.jl"))
using .Anime
using HDF5
using Plots

# We also load two instances of *observation* from existing data sets to illustrate the various models that are generated -- a single-channel data set
# and a multi-frequency data set.
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

obs1 = loadms(msname, stations, corruptseed, tropseed, tropwetonly, correff, tropattenuate, tropskynoise, tropmeandelays, tropturbulence, polvisframe,
polmode, ptginterval, ptgscale, ptgmode, gainsmode, bpfile, delim=",", ignorerepeated=false)

# Now load the multi-frequency data set

obs2 = loadms(joinpath(relativepath, "test", "data", "eht1.ms"), stations, corruptseed, tropseed, tropwetonly, correff, tropattenuate, tropskynoise, tropmeandelays, tropturbulence, polvisframe,
polmode, ptginterval, ptgscale, ptgmode, gainsmode, bpfile, delim=",", ignorerepeated=false)

# ## Atmospheric models
# At mm-wavelengths (230 GHz), the troposphere has significant effects on signal propagation. `Anime` re-implements in `Julia` all the tropospheric effects
# simulated by `MEQSv2`[^1]. The advantage here is that apart from being faster and not requiring a regular grid of complex visibilities in 
# baseline-time space, they can be called from other imaging and calibration packages when necessary.

# The function [`troposphere!`](@ref Anime.troposphere!) computes various tropospheric effects based on the flags set when `loadms` is called.

# !!! note
#     The external program `AATM` is requred to generate quantities related to tropospheric absorption and dispersion.
#     Here we pass pre-existing CSV output from `AATM` to `absorptionfile` and `dispersivefile`. The elevation angles of all antennas
#     during the course of the observation can also be preloaded when `casatools` are unavailable to generate them from scratch.

h5file = "sample.h5"
absorptionfile = joinpath(relativepath, "test", "data", "absorption1.csv")
dispersivefile = joinpath(relativepath, "test", "data", "dispersive1.csv")
elevfile = joinpath(relativepath, "test", "data", "insmodeluvf.h5")

troposphere!(obs1, h5file, absorptionfile=absorptionfile, dispersivefile=dispersivefile, elevfile=elevfile)

# This computes the delays introduced by the mean and turbulent components of the troposphere, along with attenuation due to opacity and increase in
# system noise.

# `Anime` provides plotting functions that help visualize these models. In the following, the gaps in the plotted curves signify lags between two
# observing scans.

# For example, to plot the elevation angles by station we can just do
plotelevationangle(elevfile, obs1.scanno, obs1.times, obs1.stationinfo.station, save=false)

# The transmission values computed can be plotted using
plottransmission(h5file, obs1.stationinfo.station, obs1.times, obs1.chanfreqvec, save=false)
# Since this is a channel-averaged data set, the frequency-dependent transmission reduces to a single curve per station.

# The delays due to the mean component of the troposphere can be plotted as follows:
plotmeandelays(h5file, obs1.stationinfo.station, obs1.times, obs1.chanfreqvec, save=false)
#-
rm(h5file) # hide
rm("atm.csv") # hide

# ## Primary beam response
# We model the effects of antenna pointing offsets caused by various mechanical or electronic effects. The size of the primary beam affects how much the errors
# in antenna pointing attenuate the complex visibilities. `Anime` models the pointing offsets as Gaussian processes which are useful for modelling smooth
# variations in pointing over time.

# The function [`pointing!`](@ref Anime.pointing!) is used to compute and apply pointing models to data.
pointing!(obs1, h5file=h5file)
# Note that this method is a shorthand for another method with multiple arguments that provides more fine-grained control over the input parameters.

# We now plot the pointing model generated:
plotpointingerrors(h5file, obs1.scanno, obs1.stationinfo.station, save=false)
# Mispointings of station LM (Large Millimeter Telescope, Mexico), the largest dish in the array, result in the largest attenuation of amplitude.
rm(h5file) # hide

# ## Instrumental polarization
# The feed receptors are designed to be sensitive to orthogonal polarization states in either circular or linear bases. Due to imperfections in the feed
# (either mechanical or electronic), the orthogonal measurements "leak" into the other feed, thereby giving rise to a multiplicative Jones matrix with
# small non-zero off-diagonal terms. This *feed error* or *leakage* matrix is also known as the D-Jones term. In practice, this term can vary with frequency.

# `Anime` generates smoothly varying frequency-dependent D-terms using Gaussian processes, taking a location and a scale parameter that determine the
# amount of leakage at each station. If the user requests to apply instrumental polarization to visibilities, they can be written out either in sky frame or
# antenna frame i.e., with or without parallactic angle de-rotation respectively.
inh5file = joinpath(relativepath, "test", "data", "insmodeluvf.h5")
instrumentalpolarization!(obs1, h5file=h5file, elevfile=inh5file, parangfile=inh5file)
#-
plotparallacticangle(h5file, obs1.scanno, obs1.times, obs1.stationinfo.station, save=false)
#-
plotdterms(h5file, obs1.stationinfo.station, obs1.chanfreqvec)
# There is only one frequency channel since this is a channel-averaged data set.

rm(h5file) # hide

# ## Receiver gains
# Temporal variations in complex receiver gainsfor both feeds are modelled independently. The amplitudes are modelled using a Gaussian process kernel
# (such as SE) while the phases are modelled using a Wiener process with a given location and scale parameters.

# We use an observation with only 2 scans to illustrate this better.
stationgains!(obs2, h5file=h5file)
#-
plotstationgains(h5file, obs2.scanno, obs2.times, obs2.exposure, obs2.stationinfo.station)
#-
rm(h5file) # hide

# Also there is a complex bandpass gain variation that is modelled by using representative bandpass amplitude values at certain frequencies across the bandwidth
# and interpolated for the missing frequency channels.
bandpass!(obs2, h5file=h5file)
#-
plotbandpass(h5file, obs2.stationinfo.station, obs2.chanfreqvec)
#-
rm(h5file) # hide

# ### References
# [^1]: Natarajan I. et al. MeqSilhouette v2 (2022) [MNRAS](https://academic.oup.com/mnras/article/512/1/490/6537429)