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
`Anime` is a Julia [@Bezanson2015] package for modelling the instrument response of very-long-baseline-interferometry (VLBI) arrays. It aims to model from first principles the effects of Earth's atmosphere and the electronic and mechanical properties of antennas that significantly affect VLBI observations at millimetre (mm) and sub-mm wavelengths, ensuring that the models accurately represent the nature and character of variability seen in real VLBI measurements. Such realistic instrument models are used to simulate synthetic data or correct for instrumental effects and reconstruct images of astronomical sources from observed data. It also aims to provide efficient handling of and seamless conversion routines between popular formats used to represent VLBI data and their metadata.

[^1]: https://julialang.org

# Statement of need
The Event Horizon Telescope (EHT) is a global mm-VLBI network that produced the first event-horizon-scale images of supermassive black holes at the centres of M87 [@M87PaperI] and the Milky Way [@SgrAPaperI] galaxies. New telescopes are added to the network regularly, while new locations are being considered for future stations. A physics-based approach to model atmospheric and instrumental characteristics is essential to design new instruments and understand existing ones. Such a "forward modelling" approach aids in developing new algorithms for processing data from sparse, heterogeneous VLBI arrays such as the EHT and the upcoming next-generation EHT (ngEHT).

In radio astronomy parlance, "calibration" refers to the act of removing atmospheric and instrumental effects from the data while "imaging" refers to the act of reconstructing the sky brightness distribution from the data. 
Because data calibration depends on knowledge of the sky brightness and imaging depends on knowledge of the calibration, these two procedures may be performed either iteratively or simultaneously depending on software capability. Realistic instrument models enable accurate calibration and reconstruction of sky brightness from the data. They are also used to generate synthetic data to test new algorithms.

`MEQSv2` [@Natarajan2022], a Python synthetic data generation package, introduces physically motivated signal corruptions to radio observations and is used by the end-to-end VLBI simulation and calibration pipeline `SYMBA` [@Roelofs2020] which calibrates the synthetic data using `rPICARD` [@Janssen2019]. While this introduces complex effects, it has an inflexible workflow that suffers from long runtimes when simulating long observing sessions.
`Anime` aims to provide a fast and flexible instrument modelling framework by taking advantage of the features offered by the Julia programming language, which combines the performance of languages such as C with the ease of development found in languages such as Python. Julia's automatic differentiation support also makes these instrument models inherently differentiable. 
Finally, `Comrade`, a Bayesian imaging software [@Tiede2022] that is being used increasingly in the EHT can potentially import models from `Anime` natively, enabling the development of an end-to-end framework for synthetic data generation and image reconstruction within the Julia ecosystem.

# Software components
The metadata for generating instrument models are loaded in-memory from various input data formats (Figure 1) and the output is stored as HDF5 files. Optional steps (denoted by dashed boxes) involve computing uncorrupted VLBI measurements (or "source coherency") from a given sky model using external software, applying instrument models to them, and converting between VLBI data storage formats.

![Components and control flow of a typical modelling run.](anime-components.png)

A generic instrument model may vary along station, time, frequency and polarization axes. `Anime` includes models for the troposphere which significantly affects signal propagation at mm-wavelengths (86 GHz and above). It also models the instrumental contribution to signal polarization, telescope tracking offsets, bandpass effects and receiver electronic gains, alongside the noise contributions from the atmosphere and receiver electronics. All time-variable effects are modelled using Gaussian Processes (GPs) [@GPML2006], with the hyperparameters chosen to statistically match the empirically measured temporal correlation structure of the quantity being modelled.

Synthetic data generation capabilities are built into `Anime`, with support for popular VLBI data storage formats such as UVFITS and Measurement Sets (MS)[^2]. The generated instrument models are applied to the uncorrupted data using the Radio Interferometer Measurement Equation (RIME) [@OMS2011]:

$$
\mathrm{V}_{pq} = \mathbf{\textit{G}}_p \left( \sum_{s} \mathbf{\textit{E}}_{sp}\, \mathrm{X}_{spq}\, \mathbf{\textit{E}}_{sq}^H \right) \mathbf{\textit{G}}_q^H,
$$

where $\mathrm{X}_{spq}$ is the source coherency observed towards source $s$ by the baseline formed by stations $p$ and $q$ in the absence of corrupting effects, $\mathbf{\textit{E}}_{sp}$ and $\mathbf{\textit{G}}_p$ are complex-valued matrices describing various propagation path effects and $\mathrm{V}_{pq}$ are the VLBI measurements known as "visibilities".

`Anime` can be run in modular or pipeline modes. In modular mode, the user imports `Anime` to compute instrument models by calling the relevant functions.
In pipeline mode, no user interaction is required to generate instrument models and apply them to an observation schedule. The following example computes and applies instrument models to a polarized ring-like astrophysical source observed by a sample EHT array with two polarization feeds labelled R and L.
```julia
using Anime
msfromuvfits("eht.uvfits", "eht.ms", "uvfits") # generate MS from UVFITS
run_wsclean("eht.ms", "polring", true, 1, 8191) # compute source coherency
# load data and observation parameters into Observation struct
obs = loadms("eht.ms", "eht.stations", 42, 12345, false, 0.88, false, true, true, true,
        "sky", "gp", 5.0, 100.0, "gp", "gp", "eht.bp", delim=",", ignorerepeated=false)
# compute instrument models
h5file = "models.h5"
troposphere!(obs, h5file)
instrumentalpolarization!(obs, h5file=h5file)
pointing!(obs, h5file=h5file)
stationgains!(obs, h5file=h5file)
bandpass!(obs, h5file=h5file)
thermalnoise!(obs, h5file=h5file)
postprocessms(obs, h5file=h5file) # write changes to disk
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
This work was supported by the Black Hole Initiative, which is funded by grants from the John Templeton Foundation (Grant #62286) and the Gordon and Betty Moore Foundation (Grant GBMF-8273) - although the opinions expressed in this work are those of the author(s) and do not necessarily reflect the views of these Foundations.

# References
