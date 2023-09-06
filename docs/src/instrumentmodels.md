# Instrument Models

VLBI enables the highest angular resolution achievable in astronomy, up to ~20 μas in the case of the Event Horizon Telescope (EHT) that produced the first ever images of a black hole in 2019. Since VLBI uses a sparse, heterogeneous array of radio telescopes situated around the planet, reconstructing images of observed astronomical sources is an ill-posed problem, and a deeper understanding of not only the astronomical source of interest but the instrument itself becomes crucial. Detailed explanation of the various instrument models that are simulated can be found in standard textbooks[^TMS] and other references scattered throughout this documentation.

## Radio Interferometer Measurement Equation
The Radio Interferometer Measurement Equation (RIME)[^HBS][^OMS] lies at the heart of modelling interferometric observations. A generic discrete RIME can be written as

```math
\mathbf{V}_{pq} = G_p \left( \sum_{s} E_{sp}\, \mathbf{X}_{spq}\, E_{sq}^H \right) G_q^H,
```

where the summation is carried out over all the sources $s$, and $\boldsymbol{E}_{sp}$ and $\boldsymbol{G}_p$ denote generic direction-dependent effects (DDEs) and direction-independent effects (DIEs) respectively. Each term is a $2\times2$ *Jones* matrix that describes any linear transformation acting on the incoming wave, and $H$ is the Hermitian conjugate.

Forward modelling the instrument consists of generating the Jones matrices in a physically meaningful way. `Anime` can currently model effects such as tropospheric absorption and emission and phase delays, instrumental polarization, primary beam attenuation due to mispointing, and complex bandpass and receiver gains. The additive noise terms due to the troposphere and thermal noise are also modelled.

## Gaussian Process Modelling
`Anime` uses Gaussian processes (GPs) to model time series. Gaussian processes are a generalization of the Gaussian probability distribution for functions (not just random variables). They define a probability distribution over possible functions that fit the data. A prior over the function space is encoded by the *kernel* function which takes two input data points and returns a measure of the covariance between them. Some kernel functions are manually implemented within `Anime`, with the *squared exponential (SE)* kernel (aka the *Gaussian* or the *radial basis function (RBF)* kernel) used as the default:

The SE kernel is defined as
```math
k_{SE}(x-x') = \sigma^2 e^{-\frac{(x-x')^2}{2\rho^2}}
```
The variance $\sigma^2$ is the scale factor and $\rho$ is the characteristic length scale that determines the smoothness of the function.

## The Earth's Atmosphere
Starting from a few GHz and above, the lowest layer of Earth's atmosphere, the troposphere, comes into play and is a significant contributor to signal corruptions at mm and sub-mm wavelengths. In `Anime` we characterize it as consisting of a "mean" component and an additional rapidly varying turbulent component, following the implementation in `MEQSv2`[^IN2022].

The mean troposphere introduces smoothly varying time delays that result in phase slopes with frequency. The turbulence in the troposphere introduces rapidly-varying "ad hoc" delays. The troposphere also absorbs radiation due to molecular transitions (rotational transitions of H$_2$O and O$_2$). This, along with a frequency-dependent component to the opacity of the troposphere results in an attenuation of visibility amplitudes. In thermodynamic equilibrium, the troposphere also emits radiation which increases the system temperature.

The path length due to the wet (H$_2$O) and non-wet troposphere is computed using [`AATM`](https://www.mrao.cam.ac.uk/~bn204/alma/atmomodel.html#aatm-download)[^JRP2001] and is used to generate the delays due to mean troposphere. Delays due to turbulence are simulated using Brownian random walk. The elevation and transmission-dependent sky noise due to increased system temperature is accounted for in the noise budget.

## Instrumental Polarization
While the two polarization feeds on an antenna nominally measure orthogonal polarization states of the electromagnetic wave, mechanical/electronic imperfections in the signal path both feeds to be receptive to a small fraction of the other polarization state. In addition, the mount type of any azimuthally mounted telescope causes the feeds to rotate with respect to the sky which needs to be corrected for. These corrections together constitute the instrumental contribution to the polarization state of the recorded signal.

These effects are captured using the P-Jones (parallactic angle rotation) and D-Jones (polarization leakage) matrices given by[^Hales2017]
```math
P_p = \begin{pmatrix} \mathrm{e}^{-j\chi_{_p}} && 0 \\ 0 && \mathrm{e}^{j\chi_{_p}} \end{pmatrix}
```

and
```math
D_p = \begin{pmatrix} 1 && d_{pR} \\ -d_{pL} && 1 \end{pmatrix}\, .
```

The observed visibilities are given by
```math
\mathbf{V} = (P_p^H) D_p (P_p) \mathbf{V}_{\rm sky} (P_q^H) D_q^H (P_q)
```
where
```math
P_p^H D_p P_p = \begin{pmatrix} 1 & \exp(2j\chi_p)d_p^R \\ \exp(-2j\chi_p)d_p^L & 1\end{pmatrix}\, .
```

## Primary Beams
The primary beams of participating stations determine the Field-of-View (FoV) of the observation. Several factors cause antennas to mispoint and modify its gain response, attenuating the measured visibility amplitudes, $|\mathbf{V}_{pq}|$. At mm-wavelengths even small errors in antenna pointing can cause significant attenuation.

The rms pointing error per station $\mathcal{P}_{\rm rms}$ values supplied by the user depend on site characteristics determined based on empirical measurements. The set of correlated per scan pointing offsets (${\rho_p}$) per pointing interval is obtained using GP, with $\mathcal{P}_{\rm rms}$ as the scale parameter and the scan length as the characteristic length. A Gaussian primary beam profile is then used to model the primary beamshapes of the antennas, which is an acceptable approximation at the few hundreds of $\mu$as FoV that we are interested in for mm-VLBI.
```math
E_p = \mathrm{e}^{\Bigg( -\frac{1}{2} \Bigg[\frac{\rho_p}{(\mathcal{P}_{{\rm FWHM}, p}/2\sqrt{2\ln 2})}\Bigg]^2 \Bigg)}
```
where
```math
\rho_p = \sqrt{\delta l_p^2 + \delta m_p^2}
```

## Bandpass and Receiver Gains
The generic receiver gain matrix can have both a time-varying and a frequency-varying component. The frequency-dependent component (the *B-Jones* matrix), aims to capture the variations in the amplitude and phase of the measured correlation coefficients across the bandwidth. Representative complex bandpass gains are used to spline interpolate over the entire frequency range independently for the two polarization feeds to construct the B-Jones terms.

The time-dependent gain matrix (often just referred to as *G-Jones* matrix) per station is generated using GP, with representative station-based scale parameters and the characteristic length set to the scan length to ensure smooth variation of gains over the duration of a scan. The amplitudes and phases of these terms are generated with different smoothness scales to capture the statistical properties of receiver gain evolution over time in real observations.

Both terms are modelled using diagonal 2 x 2 Jones terms, assuming that the polarization basis of the gain matrices is the same as that of the visibilities. Currently, this is the default and only option in which to represent the basis of the gain matrices.

## Noise Components
Apart from the sky noise contribution mentioned earlier, a receiver thermal noise model is also included in the noise budget. Station-based SEFD values are used to determine the per-visibility rms noise, $\sigma^{th}_{pq}$:
```math
\sigma_{pq} = \frac{1}{\eta} \sqrt{\frac{{\mathrm{SEFD}}_p {\mathrm{SEFD}}_q}{2\Delta \nu \tau}}
```
where $A_{\mathrm e}$ denotes the effective area of the telescope and $\eta$ comprises any relevant efficiency terms, such as the antenna aperture efficiency, $\eta_{\rm ap}$ and the correlator efficiency, $\eta_{\rm corr}$, $\Delta\nu$ is the bandwidth and $\tau$ the integration time. For standard 2-bit quantization, $\eta$ is set to 0.88. The SEFD itself is defined as
```math
\mathrm{SEFD} = \frac{2 k_{\mathrm B} T_{\mathrm{sys}}}{\eta_{\mathrm{ant}} A_{\mathrm e}}
```
where $\eta_{mathrm{ant}}$ is the antenna efficiency and $A_{\mathrm e}$ is the effective area of the antenna.


### References
[^TMS]: Thompson A.R., Moran J.M., Swenson Jr. G.W. Interferometry and Synthesis in Radio Astronomy (2017) [Springer](https://link.springer.com/book/10.1007/978-3-319-44431-4)
[^HBS]: Hamaker J.P., Bregman J.D., Sault R.J. Understanding radio polarimetry I (1996) [A&AS](https://articles.adsabs.harvard.edu/pdf/1996A%26AS..117..137H)
[^OMS]: Smirnov O.M. Revisiting the radio interferometry measurement equation I (2011) [A&A](https://www.aanda.org/articles/aa/pdf/2011/03/aa16082-10.pdf)
[^IN2022]: Natarajan I. et al. MeqSilhouette v2: spectrally resolved polarimetric synthetic data generation for EHT (2022) [MNRAS](https://academic.oup.com/mnras/article/512/1/490/6537429)
[^JRP2001]: Pardo J.R., et al. Atmospheric transmission at microwaves (ATM): an improved model for mm/submm applications (2001) [IEEE Xplore](https://ieeexplore.ieee.org/document/982447)
[^Hales2017]: Hales C., Calibration Errors in Interferometric Radio Polarimetry (2017) [The Astronomical Journal](https://iopscience.iop.org/article/10.3847/1538-3881/aa7aef)