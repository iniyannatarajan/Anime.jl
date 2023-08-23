# Instrument Modelling

VLBI enables the highest angular resolution achievable in astronomy, up to ~20 μas in the case of the Event Horizon Telescope (EHT) that produced the first ever images of a black hole in 2019. Since VLBI uses a sparse, heterogeneous array of radio telescopes situated around the planet, reconstructing images of observed astronomical sources is an ill-posed problem, and a deeper understanding of not only the astronomical source of interest but the instrument itself becomes crucial.

## Radio Interferometer Measurement Equation
The Radio Interferometer Measurement Equation (RIME)[^HBS][^OMS] lies at the heart of modelling interferometric observations. A generic discrete RIME can be written as

```math
\mathbf{V}_{pq} = G_p \left( \sum_{s} E_{sp}\, \mathbf{X}_{spq}\, E_{sq}^H \right) G_q^H,
```

where the summation is carried out over all the sources $s$, and $\boldsymbol{E}_{sp}$ and $\boldsymbol{G}_p$ denote generic direction-dependent effects (DDEs) and direction-independent effects (DIEs) respectively. Each term is a $2\times2$ *Jones* matrix that describes any linear transformation acting on the incoming wave, and $H$ is the Hermitian conjugate.

Forward modelling the instrument consists of generating the Jones matrices in a physically meaningful way. `Anime` can currently model effects such as tropospheric absorption and emission and phase delays, instrumental polarization, primary beam attenuation due to mispointing, and complex bandpass and receiver gains.

### Troposphere
The troposphere is the lowest layer of Earth's atmosphere and is a significant contributor to signal corruptions at mm-wavelengths. It can be characterized as consisting of a "mean" component and an additional rapidly varying turbulent component.

The mean troposphere introduces smoothly varying time delays that result in phase slopes with frequency. The path length due to the wet (H$_2$O) and non-wet troposphere is computed using [`AATM`](https://www.mrao.cam.ac.uk/~bn204/alma/atmomodel.html#aatm-download)[^JRP2001] and is used to generate the delays. The troposphere also introduces 

### References
[^HBS]: Hamaker J.P., Bregman J.D., Sault R.J. Understanding radio polarimetry I (1996) [A&AS](https://articles.adsabs.harvard.edu/pdf/1996A%26AS..117..137H)
[^OMS]: Smirnov O.M. Revisiting the radio interferometry measurement equation I (2011) [A&A](https://www.aanda.org/articles/aa/pdf/2011/03/aa16082-10.pdf)
[^JRP2001]: Pardo J.R., et al. Atmospheric transmission at microwaves (ATM): an improved model for millimeter/submillimeter applications (2001) [IEEE Xplore](https://ieeexplore.ieee.org/document/982447)