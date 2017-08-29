## ----message=FALSE-------------------------------------------------------
library(transformeR)
data("NCEP_Iberia_hus850", "NCEP_Iberia_psl", "NCEP_Iberia_ta850")

## ----fig.cap='The mean sea-level pressure field for the Iberian Peninsula. Type `help("<name_of_the_dataset>")` for further details.',fig.show='hold',message=FALSE----
plotClimatology(climatology(NCEP_Iberia_psl), backdrop.theme = "coastline",
                main = "Mean DJF SLP (1983-2002)")
plotClimatology(climatology(NCEP_Iberia_hus850), backdrop.theme = "coastline",
                main = "Mean DJF hus850 (1983-2002)")
plotClimatology(climatology(NCEP_Iberia_ta850), backdrop.theme = "coastline",
                main = "Mean DJF ta850 (1983-2002)")

## ------------------------------------------------------------------------
x <- makeMultiGrid(NCEP_Iberia_hus850, NCEP_Iberia_psl, NCEP_Iberia_ta850)

## ------------------------------------------------------------------------
data("VALUE_Iberia_tp")
y <- VALUE_Iberia_tp

## ----fig.cap='Precipitation observations used as predictand. Type `help("VALUE_Iberia_tp")` for further details.',message=FALSE----
plotClimatology(climatology(y), backdrop.theme = "countries", cex = 1.5,
                main = "Mean Winter daily precip (mm/day, 1983-2002")

## ---- echo=FALSE---------------------------------------------------------
library(downscaleR)

## ------------------------------------------------------------------------
out <- prepare_predictors(x = x,
                          y = y,
                          global.vars = NULL,
                          PCA = NULL,
                          local.predictors = NULL)
str(out)

## ------------------------------------------------------------------------
str(attributes(out))

## ------------------------------------------------------------------------
getVarNames(x)

## ---- message=FALSE------------------------------------------------------
out <- prepare_predictors(x = x,
                          y = y,
                          PCA = list(n.eofs = c(10,5,5))
)

## ------------------------------------------------------------------------
str(out)

## ------------------------------------------------------------------------
getVarNames(x)

## ------------------------------------------------------------------------
out <- prepare_predictors(x = x,
                          y = y,
                          global.vars = c("ta850", "psl"),
                          local.predictors = list(neigh.vars = "hus850",
                                                  n.neighs = 4,
                                                  neigh.fun = NULL)
)

## ------------------------------------------------------------------------
str(out)

## ------------------------------------------------------------------------
out <- prepare_predictors(x = x,
                          y = y,
                          global.vars = c("ta850", "psl"),
                          local.predictors = list(neigh.vars = "hus850",
                                                  n.neighs = 4,
                                                  neigh.fun = list(FUN = "mean"))
)

## ------------------------------------------------------------------------
str(out)

## ----message=FALSE-------------------------------------------------------
out <- prepare_predictors(x = x,
                          y = y,
                          global.vars = c("ta850", "psl"),
                          PCA = list(v.exp = c(.95, .95)),
                          local.predictors = list(neigh.vars = "hus850",
                                                  n.neighs = 4,
                                                  neigh.fun = NULL)
)

## ------------------------------------------------------------------------
str(out)

