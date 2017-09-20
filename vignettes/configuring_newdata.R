## ---- message=FALSE------------------------------------------------------
require(transformeR)
require(downscaleR)

## ------------------------------------------------------------------------
data("VALUE_Iberia_tp")
y <- VALUE_Iberia_tp

## ------------------------------------------------------------------------
data("NCEP_Iberia_hus850", "NCEP_Iberia_psl", "NCEP_Iberia_ta850")
mg <- makeMultiGrid(NCEP_Iberia_hus850, NCEP_Iberia_psl, NCEP_Iberia_ta850)
# Calibration period 1983-1994
x <- subsetGrid(mg, years = 1983:1994)

## ------------------------------------------------------------------------
# Simulation period 1995:2002
prediction.data <- subsetGrid(mg, years = 1995:2002) #%>% redim(member = TRUE)

## ------------------------------------------------------------------------
predictor <- prepare_predictors(x = x,
                                y = y,
                                global.vars = c("psl", "ta850"),
                                PCA = list(v.exp = c(.9, .95)),
                                local.predictors = list(neigh.vars = "hus850",
                                                        n.neighs = 4,
                                                        neigh.fun = list(FUN = "mean")))

## ------------------------------------------------------------------------
newdata <- localScaling(grid = prediction.data,
                        base = x,
                        ref = prediction.data,
                        time.frame = "monthly")

## ------------------------------------------------------------------------
str(newdata)

## ------------------------------------------------------------------------
test.data <- prepare_newdata(predictor = predictor,
                             newdata = newdata)
str(test.data)

## ------------------------------------------------------------------------
print(sessionInfo(package = c("transformeR", "downscaleR")))

