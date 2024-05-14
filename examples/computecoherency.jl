# # Compute coherency matrix

# !!! note
#     While `Anime` focuses on generating instrument models, this example is provided as a "passive" example for the user to run on their local
#     machines. This function cannot be used if [`WSClean`](https://wsclean.readthedocs.io/en/latest/index.html)
#     is not installed. Usage of external sofware will soon be replaced by native computation of source coherency matrices written in Julia.

# The Radio Interferometer Measurement Equation (RIME)[^HBS][^OMS] lies at the heart of modelling interferometric observations.
# A generic discrete RIME can be written as

# ```math
# \mathbf{V}_{pq} = G_p \left( \sum_{s} E_{sp}\, \mathbf{X}_{spq}\, E_{sq}^H \right) G_q^H,
# ```

# where the summation is carried out over all the sources $s$, and $\boldsymbol{E}_{sp}$ and $\boldsymbol{G}_p$ denote generic direction-dependent 
# effects (DDEs) and direction-independent effects (DIEs) respectively. Each term is a $2\times2$ *Jones* matrix that describes any linear transformation 
# acting on the incoming wave, and $H$ is the Hermitian conjugate. $X_{spq}$ is the source coherency matrix given by

# ```math
# X_{spq} = \mathrm{B} e^{-2\pi i (u_{pq}l + v_{pq}m + w_{pq}(n-1))}; \mathbf{u}_{pq} = \mathbf{u}_p - \mathbf{u}_q
# ```
# where $\mathrm{B}$ is the brightness matrix, $l$, $m$, $n$ are direction cosines, and $\mathbf{u}_{pq}$ are the baseline $uvw$ coordinates.

# While computing instrument models is the main aim of `Anime` we also provide a way to compute source coherency matrices by calling the external program
# `WSClean`:

# ```julia
# using Anime
# msname = "eht.ms" # path to existing Measurement Set
# skymodel = "grmhdpol" # path to existing sky model; samples are available in the `inputs/` directory in the source code
# polarized = true
# channelgroups = 1
# oversamplingfactor = 8191
#
# run_wsclean(msname, skymodel, polarized, channelgroups, oversamplingfactor)
# ```

# ### References
# [^HBS]: Hamaker J.P, Bregman J.D., Sault R.J. Understanding radio polarimetry (1996) [AAPS](https://articles.adsabs.harvard.edu/pdf/1996A%26AS..117..137H)
# [^OMS]: Smirnov O.M (2011) Revisiting the radio interferometry measurement equation [A&A](https://www.aanda.org/articles/aa/pdf/2011/03/aa16082-10.pdf)