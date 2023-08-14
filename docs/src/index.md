```@meta
CurrentModule = Anime
```

# Anime

**A**tmospheric a**N**d **I**nstrumental models in the **M**easurement **E**quation - `Anime` - is an instrument modelling framework for radio interferometry written in Julia, an open source, high performance language for scientific computing. `Anime` aims to support efficient handling of various data formats commonly used in VLBI, provide seamless conversion between these formats and a variety of data products as output by a Very Long Baseline Interferometry (VLBI) array, and integrate with other Julia-based software packages for VLBI data analysis.

Signals from astronomical sources are affected by various effects (e.g. atmospheric, mechanical, electronic) along the propagation path before they are recorded. Some of these effects can be modelled from first principles while for others a phenomenological approach that captures the statistical properties of the effect is more useful. Modelling these effects is an important step towards understanding both the astronomical source of interest and the capabilities and limitations of existing and planned instruments.

Within the EHT, two software packages are generally used for simulating interferometric observations of black holes: [eht-imaging](https://github.com/achael/eht-imaging/) and [SYMBA](https://bitbucket.org/M_Janssen/symba/src/master/) (which uses [MEQSv2](https://github.com/rdeane/MeqSilhouette/tree/focalpy38)). While `eht-imaging` follows an a posteriori approach to generate synthetic VLBI data whose statistical properties closely match real VLBI data, `SYMBA/MEQSv2` take a physics-based a priori approach to simulate certain propagation path effects. These packages complement each other: `eht-imaging` can introduce some data corruptions and generate data sets faster and `SYMBA/MEQSv2` can introduce more complex propagation path effects but at a higher computational cost due to its being a full-fledged pipeline with many moving parts implemented in `bash` and `python` that stitch together data products from different "monolothic" stages.

`Anime` is designed to provide a complete framework in Julia to compute instrument models from first principles wherever possible and to generate synthetic VLBI data sets with support for and conversion between multiple VLBI data storage formats. It can construct instrument models for VLBI observations with multiple scans, with irregulary-spaced and missing data without having to construct a regular grid of complex visibilities in baseline and time, thereby speeding up the generation of synthetic data sets. It also constructs more realistic models, improving upon some propagation path effects simulated by `SYMBA/MEQSv2`. The goal is to be able to read from and write to data formats that are commonly used in VLBI data analysis and provide a set of instrument modelling functions that other software written in Julia can import and use.

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