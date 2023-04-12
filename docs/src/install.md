# Installation

## Pre-requisites
Anime.jl uses the python package `casatools` to handle the creation and manipulation of [CASA Measurement Sets](https://casa.nrao.edu/Memos/229.html). If you install Anime.jl using Julia's package manager, `casatools` will be automatically installed.

Install [`WSClean`](https://wsclean.readthedocs.io/en/latest/) for predicting uncorrupted visibilities. On Ubuntu, this can be done via the Ubuntu package manager.

!!! note
    The use of `WSClean` will be deprecated soon as source coherency will be computed internally in Julia instead of relying on external packages.

Install [`AATM`](https://www.mrao.cam.ac.uk/~bn204/alma/atmomodel.html#aatm-download) for computing certain atmospheric quantities that affect the observation, such as transmission, dry and wet path lengths, and sky temperature.

!!! note
    The use of `AATM` will be deprecated soon as more advance atmospheric modelling frameworks are integrated into Anime.

## Installing Anime
Anime.jl can be installed using Julia's package manager by entering the Julia REPL and typing
```julia
using Pkg
Pkg.add("Anime")
```
or by entering package mode by typing `]` in the Julia REPL and then typing
```julia
add Anime
```