---
title: 'Anime: Atmospheric and Instrument Models in the Measurement Equation'
tags:
  - julia
  - astronomy
  - radio astronomy
  - vlbi
  - instrument modelling
authors:
  - name: Iniyan Natarajan
    orcid: 0000-0001-8242-4373
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Lindy Blackburn
    orcid: 0000-0002-9030-642X
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Paul Tiede
    orcid: 0000-0003-3826-5648
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Dominic Pesce
    orcid: 0000-0002-5278-9221
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Freek Roelofs
    orcid: 0000-0001-5461-3687
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Center for Astrophysics | Harvard & Smithsonian
   index: 1
 - name: Black Hole Initiative at Harvard University
   index: 2
date: 19 Dec 2023
bibliography: paper.bib

---

# Summary
`Anime` is a Julia [@Bezanson2015] package for modelling the instrument response of very-long-baseline-interferometry (VLBI) arrays. It aims to model from first principles the effects of Earth's atmosphere and the electronic and mechanical properties of antennas that significantly affect VLBI observations at millimetre (mm) and sub-mm wavelengths, ensuring that the models accurately represent the nature and character of variability seen in real VLBI measurements. Such realistic instrument models are used to simulate synthetic data or correct for instrumental effects during the reconstruction of images of astronomical sources from observed data. `Anime` also facilitates efficient handling of and conversion between popular formats used to represent VLBI data and their metadata.

[^1]: https://julialang.org

# Statement of need
The Event Horizon Telescope [EHT, @M87PaperII] is a global mm-VLBI network that produced the first event-horizon-scale images of supermassive black holes at the centres of M87 [@M87PaperI] and the Milky Way [@SgrAPaperI] galaxies. Existing radio telescopes are regularly added to the network, while new geographic locations are also being considered for future stations. A physics-based approach to model atmospheric and instrumental characteristics is essential to optimize the design of new instruments and to characterize existing ones. Such a "forward modelling" approach aids in developing new algorithms for processing data from sparse, heterogeneous VLBI arrays such as the EHT and future next-generation upgrades to the EHT [ngEHT, @Doeleman2023].

In radio astronomy parlance, "calibration" refers to the act of removing atmospheric and instrumental effects from the data while "imaging" refers to the act of reconstructing the sky brightness distribution from the interferometric data. While some calibration parameters may be constrained a priori, accurate calibration parameters are often fitted to the data directly through a process of "self-calibration" [@Pearson1984], which may be performed either iteratively or simultaneously with imaging depending on software capability. Realistic instrument models thus enable both statistically efficient and accurate calibration, leading to a more faithful reconstruction of sky brightness from the data.

In addition to calibration, instrument models are used in the generation of synthetic interferometric data products. `MEQSv2` [@Natarajan2022], a Python synthetic data generation package, introduces physically motivated signal corruptions to radio observations and is used by the end-to-end VLBI simulation and calibration pipeline `SYMBA` [@Roelofs2020]. While `SYMBA` is able to simulate complex effects, it has a complicated multi-stage workflow that results in long runtimes when simulating large observing sessions. `Anime` aims to provide a fast and flexible instrument modelling framework by taking advantage of the features offered by the Julia programming language, which combines the performance of languages such as C with the ease of development found in languages such as Python. Julia's automatic differentiation support also makes these instrument models inherently differentiable. 
Finally, `Comrade`, a Bayesian imaging software [@Tiede2022] that is being used increasingly in the EHT can potentially import models from `Anime` natively, enabling the development of an end-to-end framework for synthetic data generation and image reconstruction within the Julia ecosystem.

# Software components
The metadata for generating instrument models are loaded in-memory from various input data formats (Figure 1) and the output is stored as HDF5 files. Optional steps (enclosed within dashed boxes) involve computing uncorrupted VLBI measurements (or "source coherency") from a given sky model using external software, applying instrument models to them, and converting between VLBI data storage formats.

![Components and control flow of a typical modelling run.](anime-components.png)

A generic instrument model may vary along station, time, frequency and polarization axes. `Anime` includes models for the troposphere which significantly affects signal propagation at mm-wavelengths (86 GHz and above). It also models the instrumental contribution to signal polarization, telescope tracking offsets, bandpass effects and receiver electronic gains, alongside the noise contributions from the atmosphere and receiver electronics. All time-variable effects are modelled using Gaussian Processes [GPs, @GPML2006], with the hyperparameters chosen to statistically match the empirically measured temporal correlation structure of the quantity being modelled.

Synthetic data generation capabilities are built into `Anime`, with support for popular VLBI data storage formats such as UVFITS and Measurement Sets (MS)[^2]. The generated instrument models are applied to the uncorrupted data using the Radio Interferometer Measurement Equation (RIME) that describes the full polarization state of the signal with the 2x2 _Jones matrix_ formalism [@OMS2011]:

$$
\mathrm{V}_{pq} = \mathbf{\textit{G}}_p \left( \sum_{s} \mathbf{\textit{E}}_{sp}\, \mathrm{X}_{spq}\, \mathbf{\textit{E}}_{sq}^H \right) \mathbf{\textit{G}}_q^H,
$$

where $\mathrm{X}_{spq}$ is the source coherency observed towards source $s$ by the baseline formed by stations $p$ and $q$ in the absence of corrupting effects, $\mathbf{\textit{E}}_{sp}$ and $\mathbf{\textit{G}}_p$ are complex-valued matrices describing various propagation path effects and $\mathrm{V}_{pq}$ are the VLBI measurements known as "visibilities".

`Anime` can be run in modular or pipeline modes. In modular mode, the user imports `Anime` to compute instrument models by calling the relevant functions.
In pipeline mode, no user interaction is required to generate instrument models and apply them to an observation schedule. The following example computes and applies instrument models to a polarized ring-like astrophysical source observed by a sample EHT array with two polarization feeds (R and L).
```julia
# Example synthetic data generation with Anime
using Anime
createmsfromuvfits("eht.uvfits", "eht.ms", "uvfits") # generate MS from UVFITS
run_wsclean("eht.ms", "polring", true, 1, 8191) # compute source coherency
obs = readms("eht.ms") # load MS into memory
stationinfo = readstationinfo("eht.stations") # load CSV file containing station and site information
obsconfig = readobsconfig("config.yaml") # load YAML file with observation schedule information
# compute instrument models
h5file = "models.h5"
troposphere!(obs, stationinfo, obsconfig, h5file)
instrumentalpolarization!(obs, stationinfo, obsconfig, h5file=h5file)
pointing!(obs, stationinfo, obsconfig, h5file=h5file)
stationgains!(obs, stationinfo, obsconfig, h5file=h5file)
bpinfo = readbandpassinfo("eht.bandpass") # load CSV file with bandpass gain information at representative frequencies
bandpass!(obs, stationinfo, obsconfig, bpinfo, h5file=h5file)
thermalnoise!(obs, stationinfo, obsconfig, h5file=h5file)
writems(obs, h5file=h5file) # write changes to disk
```
![Amplitudes of the four polarization products with instrument models applied versus the baseline length between pairs of stations.](datavis_visampvspbs.png)

[^2]: https://casa.nrao.edu/Memos/229.html

# Related Packages
- `MEQSv2` [@Natarajan2022]: A synthetic data generation package for VLBI written in Python. It was the first VLBI simulator used in the EHT to include atmospheric effects and can compute most instrument models found in `Anime`.
- `SYMBA` [@Roelofs2020]: An end-to-end synthetic data generation pipeline that uses `MEQSv2` to generate synthetic data and introduces residual calibration effects to closely match the properties of real data.
- `eht-imaging` [@Chael2018]: A general-purpose python package for analyzing EHT observations, with simulation modules for generating synthetic data.
- `ngehtsim`: A fast and flexible (sub-)mm VLBI synthetic data generator based on `eht-imaging`, adding capabilities such as simulation of local weather effects and fringe-finding residuals.
- `Comrade` [@Tiede2022]: A Bayesian imaging framework for reconstructing images from VLBI observations while accounting for calibration residuals.

# Acknowledgements
The authors thank Dominic Chang, Alexander Plavin, and Torrance Hodgson for helpful discussions. Support for this work was provided by the NSF (AST-1935980, AST-2034306) and by the Gordon and Betty Moore Foundation through grant GBMF-10423. This work was supported by the Black Hole Initiative, which is funded by grants from the John Templeton Foundation (Grant #62286) and the Gordon and Betty Moore Foundation (Grant GBMF-8273), although the opinions expressed in this work are those of the author(s) and do not necessarily reflect the views of these Foundations.

# References
