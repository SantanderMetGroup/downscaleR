cat("***WELCOME TO downscaleR!***")
require(abind)
# Load functions
#rfuncs <- normalizePath(list.files(full.names=TRUE)[-grep("\\.md$|^init\\.R|loadGridDataset", list.files())])
exceptions <- c("ignore|\\.md$|^init\\.R|\\.Rproj|datasets|commonR|loadPredictors.GridDataset.R")
rfuncs <- list.files(full.names=TRUE)[-grep(exceptions, list.files())]
for (i in 1:length(rfuncs)) {
    message("Sourcing", rfuncs[i], sep = " ")
    source(rfuncs[i])
}
rm(list=c("rfuncs", "i", "exceptions"))
## Read vocabulary
#vocabulary <- read.csv("./init/vocabulary.csv")
#rm(list = c("i","rfuncs"))
library(utils)
if (isTRUE("sp" %in% installed.packages()[ ,1])) {
    require(sp)
} else {
    install.packages("sp", repos="http://cran.es.r-project.org", dependencies=TRUE)
}
# Init JVM
if (!isTRUE("rJava" %in% installed.packages()[ ,1])) {
    install.packages("rJava", repos = "http://cran.es.r-project.org", dependencies = TRUE)
}
require(rJava)
options(java.parameters = "-Xmx2g")
.jinit("./inst/java/netcdfAll-4.3.jar")
cat("JVM successfuly initialized \n")
cat("Ready")
#End
