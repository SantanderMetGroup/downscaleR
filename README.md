# What is downscaleR?

downscaleR is an R package is an R package for empirical-statistical downscaling focusing on daily data and covering the most popular approaches (bias correction, Model Output Statistics, Perfect Prognosis) and techniques. This package has been conceived to work in the framework of both seasonal forecasting and climate change studies. Thus, it considers ensemble members as a basic dimension of the data structure. Find out more about this package at the [downscaleR wiki](https://github.com/SantanderMetGroup/downscaleR/wiki). 

This package is part of the [climate4R bundle](http://www.meteo.unican.es/climate4r), formed by `loadeR`, `transformeR`, `downscaleR` and `visualizeR`.

The recommended installation procedure is to use the `install_github` command from the devtools R package (see the installation info in the wiki):

```r
devtools::install_github(c("SantanderMetGroup/transformeR", "SantanderMetGroup/downscaleR"))
```
**IMPORTANT:** Note that `transformeR` must be previously installed on your system.
 
**NOTE:** The utilities in `transformeR` were formerly part of `downscaleR` (up to v1.3-4). Since `downscaleR` v2.0-0, these are in `transformeR` and `downscaleR` is strictly aimed to statistical downscaling and bias correction. 

---
Reference and further information: 

Cofiño et al. (2018) The ECOMS User Data Gateway: Towards seasonal forecast data provision and research reproducibility in the era of Climate Services. **Climate Services**, http://dx.doi.org/10.1016/j.cliser.2017.07.001.
