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

### I/O
```@docs
Anime.readms
Anime.readobsconfig
Anime.readstationinfo
Anime.readbandpassinfo
Anime.readalistv5
Anime.readalistv6
Anime.createmsfromconfig
Anime.createmsfromuvfits
Anime.createuvfitsfromms
Anime.writems
Anime.makecasaantennatable
Anime.addweightcolumns
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
Anime.MeasurementSet
```

## Internal
```@docs
Anime.copymodeltodata
Anime.computeweights!
Anime.run_aatm
Anime.compute_transmission!
Anime.attenuate!
Anime.compute_skynoise!
Anime.compute_meandelays!
Anime.compute_turbulence!
Anime.squaredexponentialkernel
Anime.rationalquadratickernel
```