## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                 message = FALSE,
                 warning = FALSE,
                 fig.align = "center",
                 tidy = TRUE,
                 cache = FALSE,
                 fig.width = 9,
                 fig.height = 5,
                 fig.path = "bias_correction_figs/")

## ------------------------------------------------------------------------
# If we have not installed the "devtools" library we should do it:
# install.packages("devtools")
# library(devtools)
# Installing the downscaleR R-package
# devtools::install_github("SantanderMetGroup/downscaleR")
library(downscaleR)

## ------------------------------------------------------------------------
# precipitation
data(VALUE_Iberia_pr)
y <- VALUE_Iberia_pr
# temperature
data(VALUE_Iberia_tas)
y.t <- VALUE_Iberia_tas

## ------------------------------------------------------------------------
# precipitation
data(NCEP_Iberia_pr)
x <- gridArithmetics(NCEP_Iberia_pr, 86400, operator = "*")

# temperature
data(NCEP_Iberia_tas)
x.t <- NCEP_Iberia_tas


## ------------------------------------------------------------------------
library(visualizeR)
library(sp)
spatialPlot(climatology(x), backdrop.theme = "coastline", 
            sp.layout = list(list(SpatialPoints(getCoordinates(y)), 
                                  pch = 17, col = "black", cex = 1.5)))

## ------------------------------------------------------------------------
cal <- biasCorrection(y = y, x = x, newdata = x, precipitation = TRUE, method = "delta")

## ------------------------------------------------------------------------
cal.t <- biasCorrection (y = y.t, x = x.t, newdata = x.t, method = "delta")

## ------------------------------------------------------------------------
# precipitation
quickDiagnostics(y, x, cal, type = "daily", location = c(-2.0392, 43.3075))

## ------------------------------------------------------------------------
# temperature
quickDiagnostics(y.t, x.t, cal.t, type = "daily", location = c(-2.0392, 43.3075))

## ---- message=FALSE, warning=FALSE, results= 'hide'----------------------
# time series plotting for Igueldo station
Igueldo.cal <- subsetGrid(cal, station.id = "000234")
Ig.coords <- getCoordinates(Igueldo.cal)
Igueldo.raw <- subsetGrid(x, lonLim = Ig.coords$x, latLim = Ig.coords$y) 
temporalPlot(Igueldo.cal, Igueldo.raw, cols = c("blue", "red"))

## ---- message=FALSE, warning=FALSE, results= 'hide'----------------------
# time series plotting for Igueldo station
temporalPlot(Igueldo.cal, Igueldo.raw, cols = c("blue", "red"), x.axis = "index")

## ------------------------------------------------------------------------
#For precipitation
cal <- biasCorrection(y = y, x = x, newdata = x, precipitation = TRUE, method = "scaling", 
                       scaling.type = "multiplicative")
#For temperature
cal.t <- biasCorrection(y = y.t, x = x.t, newdata = x.t, method = "scaling",
                         scaling.type = "additive")


## ------------------------------------------------------------------------
# precipitation
quickDiagnostics(y, x, cal, type = "daily", location = c(-2.0392, 43.3075))

## ------------------------------------------------------------------------
# temperature
quickDiagnostics(y.t, x.t, cal.t, type = "daily", location = c(-2.0392, 43.3075))

## ------------------------------------------------------------------------

cal <- biasCorrection(y = y, x = x, newdata = x, precipitation = TRUE, 
                       method = "eqm")
cal.t <- biasCorrection(y = y.t, x = x.t, newdata = x.t, 
                         method = "eqm")

## ------------------------------------------------------------------------
# precipitation
quickDiagnostics(y, x, cal, type = "daily", location = c(-2.0392, 43.3075))

## ------------------------------------------------------------------------
# temperature
quickDiagnostics(y.t, x.t, cal.t, type = "daily", location = c(-2.0392, 43.3075))

## ------------------------------------------------------------------------

cal <- biasCorrection(y = y, x = x, newdata = x, precipitation = TRUE, 
                      method = "pqm", fitdistr.args = list(densfun = "gamma"))

## ------------------------------------------------------------------------
cal.t <- biasCorrection(y = y.t, x = x.t, newdata = x.t, precipitation = TRUE, 
                      method = "pqm", fitdistr.args = list(densfun = "normal"))


## ------------------------------------------------------------------------
# precipitation
quickDiagnostics(y, x, cal, type = "daily", location = c(-2.0392, 43.3075))

## ------------------------------------------------------------------------
# precipitation
quickDiagnostics(y.t, x.t, cal.t, type = "daily", location = c(-2.0392, 43.3075))

## ------------------------------------------------------------------------

cal <- biasCorrection(y = y, x = x, newdata = x, precipitation = TRUE, 
                       method = "gpqm", theta = .70)

## ------------------------------------------------------------------------
# precipitation
quickDiagnostics(y, x, cal, type = "daily", location = c(-2.0392, 43.3075))

## ------------------------------------------------------------------------
# precipitation
cal <- biasCorrection(y = y, x = x,
                            precipitation = TRUE,
                            method = "eqm",
                            window = c(30, 15),
                            wet.threshold = 0.1,
                            cross.val = "kfold",
                            folds = list(1983:1989, 1990:1996, 1997:2002))

## ------------------------------------------------------------------------
# precipitation
quickDiagnostics(y, x, cal, type = "daily", location = c(-2.0392, 43.3075))

## ------------------------------------------------------------------------
# precipitation
cal <- biasCorrection(y = y, x = x,
                            precipitation = TRUE,
                            method = "eqm",
                            window = c(30, 15),
                            wet.threshold = 0.1,
                            cross.val = "loo")

## ------------------------------------------------------------------------
# precipitation
data(VALUE_Iberia_pr)
y <- VALUE_Iberia_pr
data(NCEP_Iberia_pr)
x <- gridArithmetics(NCEP_Iberia_pr, 86400, operator = "*")

# NO window
cal <- biasCorrection(y = y, x = x,
                            precipitation = TRUE,
                            method = "eqm",
                            wet.threshold = 0.1)

# Window: calibration window of 30 days to correct each 15 day time step
cal.win <- biasCorrection(y = y, x = x,
                            precipitation = TRUE,
                            method = "eqm",
                            window = c(30, 15),
                            wet.threshold = 0.1)

## ------------------------------------------------------------------------
# NO window
quickDiagnostics(y, x, cal, type = "daily", location = c(-2.0392, 43.3075))

## ------------------------------------------------------------------------
# Window
quickDiagnostics(y, x, cal.win, type = "daily", location = c(-2.0392, 43.3075))

## ---- message=FALSE, warning=FALSE, results= 'hide'----------------------
# time series plotting for Igueldo station
Igueldo.cal <- subsetGrid(cal, station.id = "000234")
Igueldo.cal.win <- subsetGrid(cal.win, station.id = "000234")
Ig.coords <- getCoordinates(Igueldo.cal)
temporalPlot(Igueldo.cal, Igueldo.cal.win, cols = c("blue", "green"), x.axis = "index")

## ------------------------------------------------------------------------
# precipitation
data(VALUE_Iberia_pr)
y <- VALUE_Iberia_pr
data("CORDEX_Iberia_pr")
x <- CORDEX_Iberia_pr
data("CORDEX_Iberia_pr.rcp85")
newdata <- CORDEX_Iberia_pr.rcp85
cal <- biasCorrection(y = y, x = x,
                          newdata = newdata,
                          precipitation = TRUE,
                          method = "eqm",
                          extrapolation = "constant",
                          window = c(30, 15),
                          wet.threshold = 0.1)

## ---- message=FALSE, warning=FALSE, results= 'hide'----------------------
# time series plotting for Igueldo station
VALUE.obs <- subsetGrid(y, station.id = "000234") # observation
CDX.cal <- subsetGrid(cal, station.id = "000234") # CORDEX RCP8.5 corrected
CDX.hist <- interpGrid(x, getGrid(VALUE.obs)) # CORDEX historical
CDX.raw <- interpGrid(newdata, getGrid(VALUE.obs)) # CORDEX RCP8.5
temporalPlot(VALUE.obs, CDX.hist, CDX.cal, CDX.raw, 
             cols = c("black", "red", "green", "red"))

## ------------------------------------------------------------------------
# observation
data(EOBS_Iberia_tas)
y <- subsetGrid(EOBS_Iberia_tas, years = 1983:2001)

# predictor
data(CFS_Iberia_tas)
x <- subsetGrid(CFS_Iberia_tas, years = 1983:2001)

# predictor in the test period
newdata <- subsetGrid(CFS_Iberia_tas, years = 2002)

## ------------------------------------------------------------------------
library(visualizeR)
spatialPlot(climatology(y, clim.fun = list(FUN = mean, na.rm = T)), 
                backdrop.theme = "countries", 
                scales = list(draw = T))

## ------------------------------------------------------------------------
spatialPlot(climatology(x, clim.fun = list(FUN = mean, na.rm = T)), 
                backdrop.theme = "countries", 
                scales = list(draw = T))

## ------------------------------------------------------------------------
cal <- biasCorrection(y = y,
                      x = x,
                      newdata = x,
                      method = "eqm")


## ------------------------------------------------------------------------
loc <- c(-5, 42)
quickDiagnostics(y, x, cal, location = loc)

## ------------------------------------------------------------------------
spatialPlot(climatology(cal))

## ------------------------------------------------------------------------
tercilePlot(CFS_Iberia_tas, obs = EOBS_Iberia_tas, year.target = 2002, color.pal = "ypb")

## ------------------------------------------------------------------------
# In the case of the temperature these are two of the bias correction methods available:
cal1 <- biasCorrection(y = y,
                      x = x,
                      newdata = newdata, 
                      method = "scaling",
                      scaling.type = "multiplicative")
cal2 <- biasCorrection(y = y,
                      x = x,
                      newdata = newdata,
                      method = "eqm")


## ------------------------------------------------------------------------
temporalPlot(y, x, newdata, cal1, cols = c("black", "red", "red", "blue"),
             xyplot.custom = list(ylim = c(-2, 18)))

## ------------------------------------------------------------------------
temporalPlot(newdata, cal1, cols = c("red", "blue"), 
             xyplot.custom = list(ylim = c(-2, 18)))

## ------------------------------------------------------------------------
temporalPlot(cal1, cal2, cols = c("blue", "green"), 
             xyplot.custom = list(ylim = c(-2, 18)))

## ------------------------------------------------------------------------
print(sessionInfo())


## ------------------------------------------------------------------------
delta <- function(o, p, s){
      corrected <- o + (mean(s) - mean(p))
      return(corrected)
}

## ------------------------------------------------------------------------
scaling <- function(o, p, s, scaling.type){
      if (scaling.type == "additive") {
            s - mean(p) + mean(o)
      } else if (scaling.type == "multiplicative") {
            (s/mean(p)) * mean(o)
      }
}

## ------------------------------------------------------------------------
#' @title Scaling method for bias correction
#' @description Implementation of Scaling method for bias correction 
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param s A vector containing the simulated climate for the variable used in \code{x}, but considering the test period.
#' @param scaling.type Character indicating the type of the scaling method. Options are \code{"additive"} (default)
#' or \code{"multiplicative"} (see details). This argument is ignored if \code{"scaling"} is not selected as the bias correction method.
#' @author S. Herrera and M. Iturbide
scaling <- function(o, p, s, scaling.type){
      if (scaling.type == "additive") {
            s - mean(p) + mean(o)
      } else if (scaling.type == "multiplicative") {
            (s/mean(p)) * mean(o)
      }
}
#end

## ------------------------------------------------------------------------
#' @importFrom stats pgamma qgamma

## ------------------------------------------------------------------------
print(sessionInfo())

