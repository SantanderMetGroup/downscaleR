downscaleR
==========

**downscaleR** is an R package for climate analysis, with a special focus on empirical-statistical downscaling of daily data. Find out more about this package at the [downscaleR wiki](https://github.com/SantanderMetGroup/downscaleR/wiki). `downscaleR` is fully integrated with the climate data structures provided by the [`loadeR` bundle](https://github.com/SantanderMetGroup/loadeR) for data access.

***
**NOTE**: Since v2.0-0, the data manipulation tools have moved to the new dependency [`transformeR`](https://github.com/SantanderMetGroup/transformeR), and **downscaleR** is strictly aimed to statistical downscaling and bias correction.
***

The recommended installation procedure is to use the `install_github` command from the devtools R package. 

```r
devtools::install_github(c("SantanderMetGroup/transformeR", "SantanderMetGroup/downscaleR"))
```

Alternatively, you can download the source files until the most recent stable version from the [Releases Tab](https://github.com/SantanderMetGroup/downscaleR/releases)

