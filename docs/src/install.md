# Installation

`Anime` can be installed using Julia's package manager by entering the Julia REPL and typing
```julia
using Pkg
Pkg.add("Anime")
```
or by entering package mode by typing `]` in the Julia REPL and then typing
```julia
add Anime
```
The source code of `Anime` can be found [here](https://github.com/iniyannatarajan/Anime.jl).

The python packages `casatools`, `casatasks`, and `casadata` are required to handle the creation of [Measurement Sets](https://casa.nrao.edu/Memos/229.html) (MS) and conversion between uvfits and MS formats. These are automatically installed when `Anime` is installed via `Pkg`.

## External Software

Some features of `Anime` require additional software packages to be installed. These are optional and it is entirely possible to use `Anime` without them at the cost of some functionality.

[`WSClean`](https://wsclean.readthedocs.io/en/latest/)[^1] is required to compute coherency matrices from FITS images into Measurement Sets. It is not possible to do this in `Anime` without `WSClean` since `Anime` aims to generate instrument models given observation parameters and not necessarily the original uncorrupted data. This additional functionality is provided in case the user wants to create synthetic data sets from scratch and apply instrument models to them.

On Debian-based systems `WSClean` can be installed via `apt-get`:
```console
sudo apt-get install wsclean
```

[`AATM`](https://www.mrao.cam.ac.uk/~bn204/alma/atmomodel.html#aatm-download)[^2] is required to compute atmospheric quantities required for generating the atmospheric model, such as transmission, dry and wet path lengths, and sky temperature. In the absence of `AATM`, precomputed values for these quantities can still be passed as inputs in CSV format to compute the tropospheric model. Before installing `AATM` ensure that the `boost` libraries are installed. On Debian-based system this can be done via `apt-get`:
```console
sudo apt-get install libboost-program-options-dev
```
Once `boost` is installed, `AATM` can be compiled as follows:
```console
cd /path/to/aatm-source-code
./configure --prefix=/path/to/aatm-installation
make
make install
export PATH=$PATH:/path/to/install/aatm-installation/bin
```

### References
[^1]: Offringa A. R. et al. WSCLEAN: an implementation of a fast, generic wide-field imager for radio astronomy (2014) [MNRAS](https://academic.oup.com/mnras/article/444/1/606/1010067)
[^2]: Atmospheric transmission at microwaves (ATM): an improved model for millimeter/submillimeter applications [IEEE Xplore](https://ieeexplore.ieee.org/document/982447)
