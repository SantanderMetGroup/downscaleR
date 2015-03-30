downscaleR
==========

R package for bias correction and downscaling of daily climate model outputs (both seasonal forecasting and climate predictions/projections). The package allows loading and handling observations, reanalysis and global and regional model data. Find out more about downscaleR at the [downscaleR's wiki](https://github.com/SantanderMetGroup/downscaleR/wiki).

The recommended installation procedure is to use the install_github command from the devtools R package. 
```R
￼devtools::install_github(c("SantanderMetGroup/downscaleR@stable","SantanderMetGroup/downscaleR.java@stable"))
```
Note that apart from the downscaleR package, a dependent package (downscaleR.java) is installed, containing the netCDF Java API internally used by downscaleR.

