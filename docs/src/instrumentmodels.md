# Instrument Modelling

VLBI enables the highest angular resolution achievable in astronomy, up to ~20 μas in the case of the Event Horizon Telescope (EHT) that produced the first ever images of a black hole in 2019. Since VLBI uses a sparse, heterogeneous array of radio telescopes situated around the planet, reconstructing images of observed astronomical sources is an ill-posed problem, and a deeper understanding of not only the astronomical source of interest but the instrument itself becomes crucial.

## Radio Interferometer Measurement Equation
The Radio Interferometer Measurement Equation (RIME)[^HBS][^OMS] lies at the heart of modelling interferometric observations. A generic discrete RIME can be written as

```math
\mathbf{V}_{pq} = G_p \left( \sum_{s} E_{sp}\, \mathbf{X}_{spq}\, E_{sq}^H \right) G_q^H,
```

where the summation is carried out over all the sources $s$, and $\boldsymbol{E}_{sp}$ and $\boldsymbol{G}_p$ denote generic direction-dependent effects (DDEs) and direction-independent effects (DIEs) respectively. Each term is a $2\times2$ *Jones* matrix that describes any linear transformation acting on the incoming wave, and $H$ is the Hermitian conjugate.

### References
[^HBS]: Hamaker J.P, Bregman J.D., Sault R.J. (1996) [https://articles.adsabs.harvard.edu/pdf/1996A%26AS..117..137H]  
[^OMS]: Smirnov O.M (2011) [https://www.aanda.org/articles/aa/pdf/2011/03/aa16082-10.pdf]
[^Comrade]: Tiede P. Comrade: Composable Modeling of Radio Emission [JOSS](https://joss.theoj.org/papers/10.21105/joss.04457)