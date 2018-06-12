## ----warning=FALSE, results="hide", message = FALSE----------------------
library(transformeR)
data(package = "transformeR")

## ----warning=FALSE, results="hide", message = FALSE----------------------
library(downscaleR)
# Selecting predictand (y) and predictor (x)
data("VALUE_Iberia_pr","VALUE_Iberia_tas")
y <- VALUE_Iberia_tas 
data("NCEP_Iberia_hus850", "NCEP_Iberia_psl", "NCEP_Iberia_ta850")
x <- makeMultiGrid(NCEP_Iberia_hus850, NCEP_Iberia_psl, NCEP_Iberia_ta850)
x <- scaleGrid(x, type = "standardize", spatial.frame = "field") # standardizing the predictors

## ----warning=FALSE, results="hide", message = FALSE----------------------
library(visualizeR)
spatialPlot(climatology(y), backdrop.theme = "countries", colorkey=T)
spatialPlot(climatology(x), backdrop.theme = "countries")

## ----warning=FALSE, results="hide", message = FALSE----------------------
igueldo.2000 <- subsetGrid(y,station.id = "000234",years = 2000)
temporalPlot(igueldo.2000)

## ----warning=FALSE, results="hide", message = FALSE----------------------
data <- prepareData(x = x, y = y)
data <- prepareData(x = x, y = y, 
               local.predictors = list(n = 4, vars = getVarNames(x))) 
data <- prepareData(x = x, y = y,
               spatial.predictors = list(v.exp = 0.95))
data <- prepareData(x = x, y = y, 
               spatial.predictors = list(v.exp = 0.95),
               local.predictors = list(n = 4, vars = getVarNames(x)))

## ----warning=FALSE, results="hide", message = FALSE----------------------
analog <- downscale.train(data, method = "analogs", n.analogs = 1) # The analog method

## ----warning=FALSE, results="hide", message = FALSE----------------------
igueldo.2000 <- subsetGrid(y,station.id = "000234",years = 2000)
pred.2000 <- subsetGrid(analog$pred,station.id = "000234",years = 2000)
temporalPlot(igueldo.2000, pred.2000)

## ----warning=FALSE, results="hide", message = FALSE----------------------
model <- downscale.train(data, method = "GLM",family = gaussian) # Linear regression

## ----warning=FALSE, results="hide", message = FALSE----------------------
model <- downscale.train(data, method = "NN", hidden = c(10,5), output = "linear") 

## ----warning=FALSE, results="hide", message = FALSE-----------------------------------------------
newdata <- prepareNewData(x,data)
pred  <- downscale.predict(newdata, analog)

## ----warning=FALSE, results="hide", message = FALSE-----------------------------------------------
xT <- subsetGrid(x, years = 1983:1999)  # training predictors
yT <- subsetGrid(y, years = 1983:1999)   # training predictands
data <- prepareData(xT,yT)       # preparing the data 
analog <- downscale.train(data, method = "analogs", n.analogs = 1)
xt <- subsetGrid(x, years = 2000)       # test predictors
newdata <- prepareNewData(xt,data)     # preparing the new predictors
pred  <- downscale.predict(newdata, analog)  # predicting 
# visualizing the results
yt <- subsetGrid(y,years = 2000)        
temporalPlot(pred,yt)             # plotting predictions and observations

## ----warning=FALSE, results="hide", message = FALSE-----------------------------------------------
analog.cv <- downscale.cv(x = x, y = y, method = "analogs", n.analogs = 1, folds = 5,
                 spatial.predictors = list(v.exp = 0.95),
                 local.predictors = list(n = 4, vars = getVarNames(x)))

