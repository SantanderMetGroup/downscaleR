# Load functions
#rfuncs <- normalizePath(list.files(full.names=TRUE)[-grep("\\.md$|^init\\.R|loadGridDataset", list.files())])
rfuncs <- list.files("/home/juaco/workspace/gitRepo/downscaling/", full.names=TRUE)[-grep("\\.md$|^init\\.R", list.files("/home/juaco/workspace/gitRepo/downscaling/"))]
for (i in 1:length(rfuncs)) {
      source(rfuncs[i])
}
## Read vocabulary
#vocabulary <- read.csv("./init/vocabulary.csv")
#rm(list = c("i","rfuncs"))
# Init JVM
if (isTRUE("rJava" %in% installed.packages()[ ,1])) {
      require(rJava)
	  options(java.parameters = "-Xmx4g")
      .jinit("/home/juaco/workspace/gitRepo/netcdfAll-4.3.jar")
} else {
      warning("library 'rJava' not found. Please install 'rJava' and source the init.R file again.")
}
#End
