# Installation


## Pre-requisites

Anime uses the python package `casatools` to handle the creation and manipulation of [CASA Measurement Sets](https://casa.nrao.edu/Memos/229.html).
You can install `casatools` using
```bash
pip install casatools
```
At the time of writing, `casatools` does not have a release candidate for python versions higher than 3.8. If this is not your default python environment, you may want to use `pip` in a python virtual environment or a `conda` environment that runs python 3.8 (since `casatools` does not have a `conda` installation candidate either).

Install [`WSClean`](https://wsclean.readthedocs.io/en/latest/) for predicting uncorrupted visibilities. On Ubuntu, this can be done via the Ubuntu package manager.

!!! note
    The use of `WSClean` will be deprecated soon as source coherency will be computed internally in Julia instead of relying on external packages.

Install [`AATM`](https://www.mrao.cam.ac.uk/~bn204/alma/atmomodel.html#aatm-download) for computing certain atmospheric quantities that affect the observation, such as transmission, dry and wet path lengths, and sky temperature.

!!! note
    The use of `AATM` will be deprecated soon as more advance atmospheric modelling frameworks are integrated into Anime.

## Installing Anime

If using `conda`, activate the appropriate environment and set the following environment variable in the shell before installing/using Anime
```bash
export JULIA_CONDAPKG_BACKEND="Null"
```

Anime can be installed using Julia's package manager. In the Julia REPL, type
```julia
import Pkg
Pkg.add("Anime")
```