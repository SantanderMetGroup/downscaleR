downscaleR
==========

R package for bias correction and downscaling of daily climate model outputs (both seasonal forecasting and climate predictions/projections). The package allows loading and handling observations, reanalysis and global and regional model data. Find out more about downscaleR at the [downscaleR's wiki](https://github.com/SantanderMetGroup/downscaleR/wiki).

The recommended installation procedure is to use the install_github command from the devtools R package. 
```R
￼devtools::install_github(c("SantanderMetGroup/downscaleR.java@stable","SantanderMetGroup/downscaleR@stable"))
```
Note that a package dependency (`downscaleR.java`) must be first installed. It contains the netCDF-Java API internally used by `downscaleR`.

