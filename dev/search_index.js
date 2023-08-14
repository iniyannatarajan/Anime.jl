var documenterSearchIndex = {"docs":
[{"location":"api/#Anime-API","page":"Anime API","title":"Anime API","text":"","category":"section"},{"location":"api/#Contents","page":"Anime API","title":"Contents","text":"","category":"section"},{"location":"api/","page":"Anime API","title":"Anime API","text":"Pages = [\"api.md\"]","category":"page"},{"location":"api/#Index","page":"Anime API","title":"Index","text":"","category":"section"},{"location":"api/","page":"Anime API","title":"Anime API","text":"Pages = [\"api.md\"]","category":"page"},{"location":"api/#Public-API","page":"Anime API","title":"Public API","text":"","category":"section"},{"location":"api/#Storage","page":"Anime API","title":"Storage","text":"","category":"section"},{"location":"api/","page":"Anime API","title":"Anime API","text":"Anime.msfromconfig\nAnime.msfromuvfits\nAnime.mstouvfits\nAnime.loadms\nAnime.postprocessms","category":"page"},{"location":"api/#Anime.msfromconfig","page":"Anime API","title":"Anime.msfromconfig","text":"msfromconfig(msname::String, mscreationmode::String, stations::String, casaanttemplate::String, spw_centrefreq::Array{Float64, 1}, \nspw_bw::Array{Float64, 1}, spw_channels::Array{Int64, 1}, sourcedict::Dict{String, Any}, starttime::String, exposure::Float64, scans::Int64,\nscanlengths::Array{Float64, 1}, scanlag::Float64; autocorr::Bool=false, telescopename::String=\"VLBA\", feed::String=\"perfect R L\",\nshadowlimit::Float64=1e-6, elevationlimit::String=\"10deg\", stokes::String=\"RR RL LR LL\", delim::String=\",\", ignorerepeated::Bool=false)\n\nCreate Measurement Set from observation parameters. Requires casatools.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.msfromuvfits","page":"Anime API","title":"Anime.msfromuvfits","text":"msfromuvfits(uvfits::String, msname::String, mscreationmode::String, stations::String; delim::String=\",\", ignorerepeated::Bool=false)\n\nConvert UVFITS to MS. Requires casatools and casatasks.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.mstouvfits","page":"Anime API","title":"Anime.mstouvfits","text":"mstouvfits(msname::String, uvfits::String, datacolumn::String; field::String=\"\", spw::String=\"\", antenna::String=\"\",\ntimerange::String=\"\", overwrite::Bool=false)\n\nConvert MS to UVFITS. Requires casatasks.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.loadms","page":"Anime API","title":"Anime.loadms","text":"loadms(msname::String, stations::String, corruptseed::Int64, tropseed::Int64, tropwetonly::Bool, correff::Float64, tropattenuate::Bool,\ntropskynoise::Bool, tropmeandelays::Bool, tropturbulence::Bool, polframe::String, polmode::String, ptginterval::Float64, ptgscale::Float64, ptgmode::String,\nstationgainsmode::String, bandpassfile::String; delim::String=\",\", ignorerepeated::Bool=false)\n\nLoad data and metadata from MS, station and bandpass tables and return a CjlObservation object.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.postprocessms","page":"Anime API","title":"Anime.postprocessms","text":"postprocessms(obs::CjlObservation; h5file::String=\"\")\n\nAdd weight and sigma columns, reshape as needed by MS and write to disk.\n\n\n\n\n\n","category":"function"},{"location":"api/#Coherency","page":"Anime API","title":"Coherency","text":"","category":"section"},{"location":"api/","page":"Anime API","title":"Anime API","text":"Anime.run_wsclean","category":"page"},{"location":"api/#Anime.run_wsclean","page":"Anime API","title":"Anime.run_wsclean","text":"run_wsclean(msname::String, fitsdir::String, polarized::Bool, channelgroups::Int64, osfactor::Int64)\n\nCompute source coherency matrix using WSClean and populate MS.\n\n\n\n\n\n","category":"function"},{"location":"api/#Instrument-Models","page":"Anime API","title":"Instrument Models","text":"","category":"section"},{"location":"api/","page":"Anime API","title":"Anime API","text":"Anime.troposphere!\nAnime.instrumentalpolarization!\nAnime.pointing!\nAnime.stationgains!\nAnime.bandpass!\nAnime.thermalnoise!","category":"page"},{"location":"api/#Anime.troposphere!","page":"Anime API","title":"Anime.troposphere!","text":"troposphere(obs::CjlObservation, h5file::String; absorptionfile=\"\", dispersivefile=\"\", elevfile=\"\")\n\nMain function to compute various components of the tropospheric model. The actual numerical values generated are serialized in HDF5 format.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.instrumentalpolarization!","page":"Anime API","title":"Anime.instrumentalpolarization!","text":"instrumentalpolarization!(data::Array{Complex{Float32},4}, scanno::Vector{Int32}, times::Vector{Float64}, stationinfo::DataFrame, phasedir::Array{Float64,2},\npos::Array{Float64, 2}, chanfreqvec::Vector{Float64}, polframe::String, polmode::String, antenna1::Vector{Int32}, antenna2::Vector{Int32}, \nexposure::Float64, rngcorrupt::AbstractRNG; h5file::String=\"\", elevfile::String=\"\", parangfile::String=\"\")\n\nCompute frequency-varying instrumental polarization (leakage, or \"D-Jones\" terms) and apply to data. The actual numerical values are serialized as HDF5.\n\n\n\n\n\ninstrumentalpolarization!(obs::CjlObservation; h5file::String=\"\", elevfile::String=\"\", parangfile::String=\"\")\n\nShorthand for instrumental polarization function when CjlObservation struct object is available.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.pointing!","page":"Anime API","title":"Anime.pointing!","text":"pointing!(data::Array{Complex{Float32},4}, stationinfo::DataFrame, scanno::Vector{Int32}, chanfreqvec::Vector{Float64}, \nptgint::Float64, ptgscale::Float64, ptgmode::String, exposure::Float64, times::Vector{Float64}, rngcorrupt::AbstractRNG, antenna1::Vector{Int32}, \nantenna2::Vector{Int32}, numchan::Int64; α=1.0, h5file::String=\"\")\n\nCompute pointing model and apply to data. The actual numerical values are serialized in HDF5 format.\n\n\n\n\n\npointing!(obs::CjlObservation; h5file::String=\"\")\n\nShorthand for pointing model function when CjlObservation struct object is available.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.stationgains!","page":"Anime API","title":"Anime.stationgains!","text":"stationgains!(data::Array{Complex{Float32},4}, scanno::Vector{Int32}, times::Vector{Float64}, exposure::Float64, \nstationinfo::DataFrame, mode::String, rngcorrupt::AbstractRNG, antenna1::Vector{Int32}, antenna2::Vector{Int32}, numchan::Int64; h5file::String=\"\")\n\nCompute time-variable complex station gains and apply to data. The actual numerical values are serialized in HDF5 format.\n\n\n\n\n\nstationgains!(obs::CjlObservation; h5file::String=\"\")\n\nShorthand for station gains function when CjlObservation struct object is available.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.bandpass!","page":"Anime API","title":"Anime.bandpass!","text":"bandpass!(data::Array{Complex{Float32},4}, bandpassfile::String, stationinfo::DataFrame, rngcorrupt::AbstractRNG,\nantenna1::Vector{Int32}, antenna2::Vector{Int32}, numchan::Int64, chanfreqvec::Vector{Float64}; h5file::String=\"\")\n\nCompute the complex bandpass model and apply to data. The actual numerical values are serialized in HDF5 format.\n\n\n\n\n\nbandpass!(obs::CjlObservation; h5file::String=\"\")\n\nShorthand for bandpass function when CjlObservation struct object is available.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.thermalnoise!","page":"Anime API","title":"Anime.thermalnoise!","text":"thermalnoise!(data::Array{Complex{Float32},4}, times::Vector{Float64}, antenna1::Vector{Int32}, antenna2::Vector{Int32}, correff::Float64,\nexposure::Float64, chanwidth::Float64, rngcorrupt::AbstractRNG, sefd::Vector{Float64}; h5file::String=\"\", noisefile::String=\"\")\n\nCompute per-baseline thermal noise using radiometer equation and apply to data. The actual numerical values are serialized in HDF5 format.\n\n\n\n\n\nthermalnoise!(obs::CjlObservation; h5file::String=\"\", noisefile::String=\"\")\n\nShorthand for thermal noise function when CjlObservation struct object is available.\n\n\n\n\n\n","category":"function"},{"location":"api/#Statistics","page":"Anime API","title":"Statistics","text":"","category":"section"},{"location":"api/","page":"Anime API","title":"Anime API","text":"Anime.genseries1d!","category":"page"},{"location":"api/#Anime.genseries1d!","page":"Anime API","title":"Anime.genseries1d!","text":"genseries1d!(series::Vector{ComplexF32}, mode::String, location::ComplexF32, scale::Float32, driftrate::Float32, nsamples::Int64, rng::AbstractRNG)\n\nGenerate a complex-valued Gaussian process 1-D series of length nsamples with the given location, scale, and driftrate parameters.\n\n\n\n\n\ngenseries1d!(series::Vector{Float32}, mode::String, location::Float32, scale::Float32, driftrate::Float32, nsamples::Int64, rng::AbstractRNG)\n\nGenerate a Float32-valued Gaussian process 1-D series of length nsamples with the given location, scale, and driftrate parameters.\n\n\n\n\n\ngenseries1d!(series::Vector{Float64}, mode::String, location::Float64, scale::Float64, driftrate::Float64, nsamples::Int64, rng::AbstractRNG)\n\nGenerate a Float64-valued Gaussian process 1-D series of length nsamples with the given location, scale, and driftrate parameters.\n\n\n\n\n\ngenseries1d!(series, times, rng::AbstractRNG; μ=0.0, σ=1.0, ρ=1.0)\n\nGenerate 1-D series using SE kernel.\n\n\n\n\n\ngenseries1d!(series, times, rng::AbstractRNG, α::Float64; μ=0.0, σ=1.0, ρ=2.0)\n\nGenerate 1-D series using RQ kernel.\n\n\n\n\n\n","category":"function"},{"location":"api/#Utils","page":"Anime API","title":"Utils","text":"","category":"section"},{"location":"api/","page":"Anime API","title":"Anime API","text":"Anime.elevationangle\nAnime.parallacticangle","category":"page"},{"location":"api/#Anime.elevationangle","page":"Anime API","title":"Anime.elevationangle","text":"elevationangle(times::Vector{Float64}, phasedir::Array{Float64,2}, stationinfo::DataFrame, pos::Array{Float64, 2})\n\nCompute elevation angle for all stations for all times.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.parallacticangle","page":"Anime API","title":"Anime.parallacticangle","text":"parallacticangle(times::Vector{Float64}, phasedir::Array{Float64,2}, stationinfo::DataFrame, pos::Array{Float64,2})\n\nCompute parallactic angle for all stations for all times.\n\n\n\n\n\n","category":"function"},{"location":"api/#Plotting","page":"Anime API","title":"Plotting","text":"","category":"section"},{"location":"api/","page":"Anime API","title":"Anime API","text":"Anime.plotvis\nAnime.plotuvcov\nAnime.plotstationgains\nAnime.plotbandpass\nAnime.plotpointingerrors\nAnime.plotelevationangle\nAnime.plotparallacticangle\nAnime.plotdterms\nAnime.plottransmission\nAnime.plotmeandelays","category":"page"},{"location":"api/#Anime.plotvis","page":"Anime API","title":"Anime.plotvis","text":"plotvis(uvw::Matrix{Float64}, chanfreqvec::Array{Float64,1}, flag::Array{Bool,4}, data::Array{Complex{Float32},4},\nnumchan::Int64, times::Vector{Float64}; saveprefix=\"data_\")\n\nPlot complex visibilities against time and projected baseline length.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.plotuvcov","page":"Anime API","title":"Anime.plotuvcov","text":"plotuvcov(uvw::Matrix{Float64}, flagrow::Vector{Bool}, chanfreqvec::Vector{Float64}; saveprefix=\"test_\")\n\nPlot uv-coverage of observation.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.plotstationgains","page":"Anime API","title":"Anime.plotstationgains","text":"plotstationgains(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, exposure::Float64, stationnames::Vector{String3})\n\nPlot complex station gains against time.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.plotbandpass","page":"Anime API","title":"Anime.plotbandpass","text":"plotbandpass(h5file::String, stationnames::Vector{String3}, chanfreqvec::Vector{Float64})\n\nPlot bandpass gains against time.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.plotpointingerrors","page":"Anime API","title":"Anime.plotpointingerrors","text":"plotpointingerrors(h5file::String, scanno::Vector{Int32}, stationnames::Vector{String3})\n\nPlot pointing errors.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.plotelevationangle","page":"Anime API","title":"Anime.plotelevationangle","text":"plotelevationangle(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, stationnames::Vector{String3})\n\nPlot evolution of station elevation angles during the course of the observation.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.plotparallacticangle","page":"Anime API","title":"Anime.plotparallacticangle","text":"plotparallacticangle(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, stationnames::Vector{String3})\n\nPlot evolution of station parallactic angles during the course of the observation.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.plotdterms","page":"Anime API","title":"Anime.plotdterms","text":"plotdterms(h5file::String, stationnames::Vector{String3})\n\nPlot frequency-dependent complex instrumental polarization.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.plottransmission","page":"Anime API","title":"Anime.plottransmission","text":"plottransmission(h5file::String, stationnames::Vector{String3}, times::Vector{Float64}, chanfreqvec::Vector{Float64})\n\nPlot tropospheric transmission variation with frequency.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.plotmeandelays","page":"Anime API","title":"Anime.plotmeandelays","text":"plotmeandelays(h5file::String, stationnames::Vector{String3}, times::Vector{Float64}, chanfreqvec::Vector{Float64})\n\nPlot mean delays against time.\n\n\n\n\n\n","category":"function"},{"location":"api/#Data-Types","page":"Anime API","title":"Data Types","text":"","category":"section"},{"location":"api/","page":"Anime API","title":"Anime API","text":"Anime.CjlObservation","category":"page"},{"location":"api/#Anime.CjlObservation","page":"Anime API","title":"Anime.CjlObservation","text":"struct CjlObservation{T} <: Anime.AbstractObservation{T}\n\nContainer type for storing observation parameters and data.\n\nFields\n\nmsname: Name of the Measurement Set\n\ndata: Complex visibility data (MS DATA column)\n\nflag: Flag array of the same dimensions as data (MS FLAG column)\n\nflagrow: Flag row Boolean vector (MS FLAG_ROW column)\n\nantenna1: Antenna 1 in a baseline pair (MS ANTENNA1 column)\n\nantenna2: Antenna 2 in a baseline pair (MS ANTENNA2 column)\n\nuvw: uvw coordinates (MS UVW column)\n\ntimes: Timestamps (MS TIME column)\n\nexposure: Integration time\n\nscanno: Scan numbers (MS SCAN_NUMBER column)\n\nnumchan: Number of frequency channels\n\nchanfreqvec: Channel frequencies (Hz) (MS CHAN_FREQ column)\n\nchanwidth: Width of frequency channel\n\nphasedir: Direction of phase centre\n\npos: Antenna positions in x,y,z coordinates\n\nstationinfo: Dataframe of station information input by user\n\ntropwetonly: Consider only water vapour in the troposphere\n\ncorreff: Correlator efficiency\n\ntropattenuate: Introduce attenuation by the troposphere\n\ntropskynoise: Introduce noise due to troposphere\n\ntropmeandelays: Introduce mean delays due to troposphere\n\ntropturbulence: Introduce turbulence in troposphere\n\npolframe: Polarization frame\n\npolmode: Mode to use to create frequency-varying D-term samples\n\nptginterval: Pointing interval (s)\n\nptgscale: Pointing scale mixture parameter\n\nptgmode: Mode to use to create pointing error time samples\n\nstationgainsmode: Mode to use to create station gain time samples\n\nbandpassfile: Input bandpass file\n\nrngcorrupt: RNG for all models except troposphere\n\nrngtrop: RNG for tropospheric models\n\n\n\n\n\n","category":"type"},{"location":"api/#Internal","page":"Anime API","title":"Internal","text":"","category":"section"},{"location":"api/","page":"Anime API","title":"Anime API","text":"Anime.addweightcols\nAnime.makecasaanttable\nAnime.copymodeltodata\nAnime.computeweights!\nAnime.run_aatm\nAnime.compute_transmission!\nAnime.attenuate!\nAnime.compute_skynoise!\nAnime.compute_meandelays!\nAnime.compute_turbulence!\nAnime.squaredexponentialkernel","category":"page"},{"location":"api/#Anime.addweightcols","page":"Anime API","title":"Anime.addweightcols","text":"addweightcols(msname::String, mode::String, sigmaspec::Bool, weightspec::Bool)\n\nAdd WEIGHT_SPECTRUM and SIGMA_SPECTRUM columns to MS. Requires casatools.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.makecasaanttable","page":"Anime API","title":"Anime.makecasaanttable","text":"makecasaanttable(stations::String, casaanttemplate::String; delim::String=\",\", ignorerepeated::Bool=false)\n\nCreate CASA antenna table from station info file. Requires casatools.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.copymodeltodata","page":"Anime API","title":"Anime.copymodeltodata","text":"copymodeltodata(msname::String)\n\nCopy MODEL_DATA to DATA in MS.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.computeweights!","page":"Anime API","title":"Anime.computeweights!","text":"computeweights!(totalrmsspec::Array{Float32, 4}, totalwtspec::Array{Float32, 4}; h5file::String=\"\")\n\nCompute total rms (sigma) values and inverse-squared visibility weights from thermal+sky noise terms stored in input HDF5 file.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.run_aatm","page":"Anime API","title":"Anime.run_aatm","text":"run_aatm(obs::CjlObservation; absorptionfile::String=\"\", dispersivefile::String=\"\")::DataFrame\n\nRun AATM (Bjona Nikolic; Pardo et al. 2001) to compute absorption by and dispersive delay in the troposphere. If AATM is not installed, this function can still accept input absorption and dispersion values in a specific CSV format and populate atm.csv.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.compute_transmission!","page":"Anime API","title":"Anime.compute_transmission!","text":"compute_transmission!(transmission::Array{Float64, 3}, obs::CjlObservation, atmdf::DataFrame, elevationmatrix::Array{Float64, 2}, g::HDF5.Group)::Array{Float64, 3}\n\nCompute elevation-dependent (mean) tropospheric transmission given opacity τ and elevation angle θ for each station.\n\ne^-τsintheta_rm el\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.attenuate!","page":"Anime API","title":"Anime.attenuate!","text":"attenuate!(obs::CjlObservation, transmission::Array{Float64, 3})\n\nAttenuate the signal as it passes through (mean) troposphere using precomputed transmission values.\n\nI = I_0 e^-τsintheta_rm el\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.compute_skynoise!","page":"Anime API","title":"Anime.compute_skynoise!","text":"compute_skynoise!(obs::CjlObservation, atmdf::DataFrame, transmission::Array{Float64, 3}, g::HDF5.Group)\n\nCompute sky contribution to visibility noise using the radiometer equation.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.compute_meandelays!","page":"Anime API","title":"Anime.compute_meandelays!","text":"compute_meandelays!(obs::CjlObservation, atmdf::DataFrame, elevationmatrix::Array{Float64, 2}, g::HDF5.Group)\n\nCompute delays due to (mean) troposphere.\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.compute_turbulence!","page":"Anime API","title":"Anime.compute_turbulence!","text":"compute_turbulence!(obs::CjlObservation, atmdf::DataFrame, elevationmatrix::Array{Float64, 2}, g::HDF5.Group)\n\nCompute phase delays due to tropospheric turbulence. The time series of phase errors for station p is given by\n\ndelta phi_p(t nu) = frac1sqrtsin(theta_mathrmel(t)) deltaphi^prime_p(t) big(fracnunu_0big)\n\nwhere ν is the list of channel frequencies and ν_0 is the reference frequency (lowest in the band).\n\n\n\n\n\n","category":"function"},{"location":"api/#Anime.squaredexponentialkernel","page":"Anime API","title":"Anime.squaredexponentialkernel","text":"squaredexponentialkernel(x1, x2; σ=1.0, ρ=1.0)\n\nGenerate squared exponential kernel function of the form\n\nk_SE(x-x) = sigma^2 e^-frac(x-x)^22rho^2\n\nwhere σ^2 is the variance and ρ is the characteristic length.\n\n\n\n\n\n","category":"function"},{"location":"examples/#Example-Anime-Workflow","page":"Example Anime Workflow","title":"Example Anime Workflow","text":"","category":"section"},{"location":"examples/","page":"Example Anime Workflow","title":"Example Anime Workflow","text":"Anime takes pre-generated source models and some basic information about weather conditions and the dishes on-site at each station as inputs.","category":"page"},{"location":"insmodelling/#Instrument-Modelling","page":"Instrument Modelling","title":"Instrument Modelling","text":"","category":"section"},{"location":"insmodelling/","page":"Instrument Modelling","title":"Instrument Modelling","text":"VLBI enables the highest angular resolution achievable in astronomy, up to ~20 μas in the case of the Event Horizon Telescope (EHT) that produced the first ever images of a black hole in 2019. Since VLBI uses a sparse, heterogeneous array of radio telescopes situated around the planet, reconstructing images of observed astronomical sources is an ill-posed problem, and a deeper understanding of not only the astronomical source of interest but the instrument itself becomes crucial.","category":"page"},{"location":"insmodelling/#Radio-Interferometer-Measurement-Equation","page":"Instrument Modelling","title":"Radio Interferometer Measurement Equation","text":"","category":"section"},{"location":"insmodelling/","page":"Instrument Modelling","title":"Instrument Modelling","text":"The Radio Interferometer Measurement Equation (RIME)[HBS][OMS] lies at the heart of modelling interferometric observations. A generic discrete RIME can be written as","category":"page"},{"location":"insmodelling/","page":"Instrument Modelling","title":"Instrument Modelling","text":"mathbfV_pq = G_p left( sum_s E_sp mathbfX_spq E_sq^H right) G_q^H","category":"page"},{"location":"insmodelling/","page":"Instrument Modelling","title":"Instrument Modelling","text":"where the summation is carried out over all the sources s, and boldsymbolE_sp and boldsymbolG_p denote generic direction-dependent effects (DDEs) and direction-independent effects (DIEs) respectively. Each term is a 2times2 Jones matrix that describes any linear transformation acting on the incoming wave, and H is the Hermitian conjugate.","category":"page"},{"location":"insmodelling/#References","page":"Instrument Modelling","title":"References","text":"","category":"section"},{"location":"insmodelling/","page":"Instrument Modelling","title":"Instrument Modelling","text":"[HBS]: Hamaker J.P, Bregman J.D., Sault R.J. (1996) [https://articles.adsabs.harvard.edu/pdf/1996A%26AS..117..137H]  ","category":"page"},{"location":"insmodelling/","page":"Instrument Modelling","title":"Instrument Modelling","text":"[OMS]: Smirnov O.M (2011) [https://www.aanda.org/articles/aa/pdf/2011/03/aa16082-10.pdf]","category":"page"},{"location":"insmodelling/","page":"Instrument Modelling","title":"Instrument Modelling","text":"[Comrade]: Tiede P. Comrade: Composable Modeling of Radio Emission JOSS","category":"page"},{"location":"install/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"install/","page":"Installation","title":"Installation","text":"Anime can be installed using Julia's package manager by entering the Julia REPL and typing","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"using Pkg\nPkg.add(\"Anime\")","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"or by entering package mode by typing ] in the Julia REPL and then typing","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"add Anime","category":"page"},{"location":"install/#External-Software","page":"Installation","title":"External Software","text":"","category":"section"},{"location":"install/","page":"Installation","title":"Installation","text":"Some features of Anime require external software to be installed. These are optional and it is entirely possible to use Anime without them at the cost of some functionality.","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"The python packages casatools, casatasks, and casadata are required to handle the creation of Measurement Sets (MS) and conversion between uvfits and MS formats. These are automatically installed when Anime is installed via Pkg.","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"WSClean is required to compute coherency matrices from FITS images into Measurement Sets. This is not possible to do in Anime without WSClean currently. Complex visibilities in uvfits/MS formats with precomputed coherency matrices can still be used. On Debian-based systems WSClean can be installed via apt-get:","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"sudo apt-get install wsclean","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"AATM is required for computing atmospheric quantities required for generating the atmospheric model, such as transmission, dry and wet path lengths, and sky temperature. In the absence of AATM, precomputed values for these quantities can still be passed in CSV format to compute tropospheric model. Before installing AATM ensure that the boost libraries are installed. On Debian-based system this can be done via apt-get:","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"sudo apt-get install libboost-program-options-dev","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"Once boost is installed, AATM can be compiled as follows:","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"cd /path/to/aatm-source-code\n./configure --prefix=/path/to/aatm-installation\nmake\nmake install\nexport PATH=$PATH:/path/to/install/aatm-installation/bin","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Anime","category":"page"},{"location":"#Anime","page":"Home","title":"Anime","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Atmospheric aNd Instrumental models in the Measurement Equation - Anime - is an instrument modelling framework for radio interferometry written in Julia, an open source, high performance language for scientific computing. Anime aims to support efficient handling of various data formats commonly used in VLBI, provide seamless conversion between these formats and a variety of data products as output by a Very Long Baseline Interferometry (VLBI) array, and integrate with other Julia-based software packages for VLBI data analysis.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Signals from astronomical sources are affected by various effects (e.g. atmospheric, mechanical, electronic) along the propagation path before they are recorded. Some of these effects can be modelled from first principles while for others a phenomenological approach that captures the statistical properties of the effect is more useful. Modelling these effects is an important step towards understanding both the astronomical source of interest and the capabilities and limitations of existing and planned instruments.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Within the EHT, two software packages are generally used for simulating interferometric observations of black holes: eht-imaging and SYMBA (which uses MEQSv2). While eht-imaging follows an a posteriori approach to generate synthetic VLBI data whose statistical properties closely match real VLBI data, SYMBA/MEQSv2 take a physics-based a priori approach to simulate certain propagation path effects. These packages complement each other: eht-imaging can introduce some data corruptions and generate data sets faster and SYMBA/MEQSv2 can introduce more complex propagation path effects but at a higher computational cost due to its being a full-fledged pipeline with many moving parts implemented in bash and python that stitch together data products from different \"monolothic\" stages.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Anime is designed to provide a complete framework in Julia to compute instrument models from first principles wherever possible and to generate synthetic VLBI data sets with support for and conversion between multiple VLBI data storage formats. It can construct instrument models for VLBI observations with multiple scans, with irregulary-spaced and missing data without having to construct a regular grid of complex visibilities in baseline and time, thereby speeding up the generation of synthetic data sets. It also constructs more realistic models, improving upon some propagation path effects simulated by SYMBA/MEQSv2. The goal is to be able to read from and write to data formats that are commonly used in VLBI data analysis and provide a set of instrument modelling functions that other software written in Julia can import and use.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"index.md\",\n    \"install.md\",\n    \"insmodelling.md\",\n    \"examples.md\",\n    \"api.md\"\n]","category":"page"}]
}
