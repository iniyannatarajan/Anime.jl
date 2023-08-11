# Anime API

## Contents
```@contents
Pages = ["api.md"]
```

## Index
```@index
Pages = ["api.md"]
```

## Public API

### Storage
```@docs
Anime.msfromconfig
Anime.msfromuvfits
Anime.mstouvfits
Anime.loadms
Anime.postprocessms
```

### Coherency
```@docs
Anime.run_wsclean
```

### Instrument Models
```@docs
Anime.troposphere!
Anime.instrumentalpolarization!
Anime.pointing!
Anime.stationgains!
Anime.bandpass!
Anime.thermalnoise!
```

### Statistics
```@docs
Anime.genseries1d!
```

### Utils
```@docs
Anime.elevationangle
Anime.parallacticangle
```

### Plotting
```@docs
Anime.plotvis
Anime.plotuvcov
Anime.plotstationgains
Anime.plotbandpass
Anime.plotpointingerrors
Anime.plotelevationangle
Anime.plotparallacticangle
Anime.plotdterms
Anime.plottransmission
Anime.plotmeandelays
```

## Data Types
```@docs
Anime.CjlObservation
```

## Internal
```@docs
Anime.addweightcols
Anime.makecasaanttable
Anime.copymodeltodata
Anime.computeweights!
Anime.run_aatm
Anime.compute_transmission!
Anime.attenuate!
Anime.compute_skynoise!
Anime.compute_meandelays!
Anime.compute_turbulence!
Anime.squaredexponentialkernel
```