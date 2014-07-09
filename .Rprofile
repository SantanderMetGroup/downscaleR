cat("***WELCOME TO downscaleR!***")
require(abind)
require(downscaleR.java)
# Load functions
#rfuncs <- normalizePath(list.files(full.names=TRUE)[-grep("\\.md$|^init\\.R|loadGridDataset", list.files())])
exceptions <- c("loadPredictors.GridDataset.R")
rfuncs <- list.files("R", full.names=TRUE)
rfuncs <- rfuncs[-grep(exceptions, rfuncs)]
for (i in 1:length(rfuncs)) {
    message("Sourcing ", rfuncs[i])
    source(rfuncs[i])
}
rm(list=c("rfuncs", "i", "exceptions"))


library(fields, verbose=FALSE)
plotMeanField <- function (seasonForecastObj) {
    sf <- seasonForecastObj
    dimNames <- attr(sf$Data, "dimensions")
    mar <- grep("^lon.*|^lat.*", dimNames)
    if (length(mar) != 2) {
        stop("Not a rectangular spatial domain")
    }
    aux <- apply(sf$Data, FUN = mean, MARGIN = mar)
    image.plot(sf$xyCoords$x, sf$xyCoords$y, aux, xlab = "X", ylab = "y", asp = 1)
    world(add = TRUE)
} 
# End

## Read vocabulary
#vocabulary <- read.csv("./init/vocabulary.csv")
#rm(list = c("i","rfuncs"))
# library(utils)
# if (isTRUE("sp" %in% installed.packages()[ ,1])) {
#     require(sp)
# } else {
#     install.packages("sp", repos="http://cran.es.r-project.org", dependencies=TRUE)
# }
# # Init JVM
# if (!isTRUE("rJava" %in% installed.packages()[ ,1])) {
#     install.packages("rJava", repos = "http://cran.es.r-project.org", dependencies = TRUE)
# }
# require(rJava)
# options(java.parameters = "-Xmx2g")
# .jinit("./inst/java/netcdfAll-4.3.jar")
# cat("JVM successfuly initialized \n")
# cat("Ready")
#End
# library(maptools)
# wpolys <- readShapePoly("/home/DATA//CARTOGRAPHY//BASICA//World/LandPolys/ne_110m_land")
# wlines <- as(wpolys, "SpatialLines")
