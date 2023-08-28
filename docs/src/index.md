```@meta
CurrentModule = Anime
```

# Anime

**A**tmospheric a**N**d **I**nstrumental models in the **M**easurement **E**quation - `Anime` - is an instrument modelling framework for radio interferometry written in Julia, an open source, high performance language for scientific computing. `Anime` aims to support efficient handling of various data formats commonly used in VLBI, provide seamless conversion between these formats and a variety of data products as output by a Very Long Baseline Interferometry (VLBI) array, and integrate with other Julia-based software packages for VLBI data analysis.

Signals from astronomical sources are affected by various effects (e.g. atmospheric, mechanical, electronic) along the propagation path before they are recorded. Some of these effects can be modelled from first principles while for others a phenomenological approach that captures the statistical properties of the effect is more useful. Modelling these effects is an important step towards understanding both the astronomical source of interest and the capabilities and limitations of existing and planned instruments.

Within the EHT, two software packages are generally used for simulating interferometric observations of black holes: `ngehtsim`[^DP] and `SYMBA`[^RJ2020]. While `ngehtsim` follows an a posteriori approach to generate synthetic VLBI data whose statistical properties closely match real VLBI data by adding new functionality to the core `eht-imaging`[^AC2018] functions, `SYMBA` uses `MeqSilhouette v2 (MEQSv2)`[^IN2022][^TB2017] to generate uncalibrated synthetic data using a physics-based a priori approach to simulate certain propagation path effects and `rPICARD`[^MJ2019] to calibrate these data to simulate realistic calibration errors. These packages complement each other: `ngehtsim` can introduce some data corruptions and generate data sets faster and `SYMBA/MEQSv2` can introduce more complex propagation path effects but at a higher computational cost due to its being a full-fledged pipeline with many moving parts implemented in `bash` and `python` that stitch together data products from different "monolothic" stages.

`Anime` is designed to provide a complete framework in Julia to compute instrument models from first principles wherever possible and to generate synthetic VLBI data sets with support for and convert b/w multiple VLBI data storage formats. It can construct instrument models for multi-scan VLBI observations with irregulary-spaced and missing data without having to construct a regular grid of complex visibilities in baseline and time (like `SYMBA/MEQSv2`), thereby speeding up the generation of synthetic data. It also constructs more realistic models, re-implementing existing models from `SYMBA/MEQSv2` in `Julia` for more efficient computation and introducing new ones. The ultimate aim is to be able to read/write commonly used on-disk data formats in VLBI and provide a set of instrument modelling functions that other software written in `Julia` (e.g. Comrade[^PT2022]) can import and use for synthetic data generation, calibration and imaging. Inference methods for emulating the statistical properties of actual data are also in development.

```@contents
Pages = [
    "index.md",
    "install.md",
    "components.md",
    "instrumentmodels.md",
    "examples.md",
    "api.md"
]
```

### References
[^DP]: Pesce et al. (in prep) [GitHub](https://github.com/Smithsonian/ngehtsim) ([Docs](https://smithsonian.github.io/ngehtsim/html/docs/source/index.html))
[^RJ2020]: Roelofs F., Janssen M., Natarajan I. et al. SYMBA: An end-to-end VLBI synthetic data generation pipeline (2020) [A&A](https://www.aanda.org/articles/aa/full_html/2020/04/aa36622-19/aa36622-19.html)
[^AC2018]: Chael A. et al., Interferometric Imaging Directly with Closure Phases and Closure Amplitudes (2018) [ApJ](https://iopscience.iop.org/article/10.3847/1538-4357/aab6a8)
[^IN2022]: Natarajan I. et al. MeqSilhouette v2: spectrally resolved polarimetric synthetic data generation for EHT (2022) [MNRAS](https://academic.oup.com/mnras/article/512/1/490/6537429)
[^TB2017]: Blecher T. et al. MEQSILHOUETTE: a mm-VLBI observation and signal corruption simulator (2017) [MNRAS](https://academic.oup.com/mnras/article/464/1/143/2194682)
[^MJ2019]: Janssen M. et al. rPICARD: A CASA-based calibration pipeline for VLBI data (2019) [A&A](https://www.aanda.org/articles/aa/full_html/2019/06/aa35181-19/aa35181-19.html)
[^PT2022]: Tiede P. Comrade: Composable Modeling of Radio Emission (2022) [JOSS](https://joss.theoj.org/papers/10.21105/joss.04457)