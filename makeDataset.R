# TODO: include argument pattern for inclusion of selected files according to regular expressions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


makeDataset <- function(source.dir, ncml.file, pattern = NULL, recursive = FALSE) {
      source.dir <- source.dir
      ncml.file <- ncml.file
      pattern <- pattern
      recursive <- recursive
      message(paste("[",Sys.time(),"]", sep=""))
      file.create(ncml.file)
      if (is.null(pattern)) {
            lf <- normalizePath(list.files(source.dir, full.names = TRUE, pattern = "\\.nc$", recursive = recursive), winslash = "/")      
      } else {
            lf <- normalizePath(list.files(source.dir, full.names = TRUE, pattern = pattern, recursive = recursive), winslash = "/")      
            lf <- lf[grep("\\.nc$", lf)]
      }
# 	lf <- normalizePath(list.files(source.dir, full.names = TRUE, pattern = pattern, recursive = recursive), winslash = "/")
	message(paste("Creating dataset from", length(lf), "files"))
	z <- file(ncml.file, "w")
      # [http://www.unidata.ucar.edu/software/netcdf/ncml/v2.2/Aggregation.html]
      cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file = z)
      cat(c("\n","<netcdf xmlns=\"http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2\">"), sep = "", file = z)
      cat(c("\n","\t","<aggregation type=\"union\">"), sep = "", file = z)
      varNames <- rep(NA, length(lf))
      for (i in 1:length(lf)) {
            jString <- .jnew("java/lang/String", lf[i])
            gds <- J("ucar.nc2.dt.grid.GridDataset")$open(jString)
            varNames[i] <- gds$getGrids()$toString()
            gds$close()
      }
      vars <- unique(varNames)
	      for (i in (1:length(vars))) {
            subset <- paste("<netcdf location=", "\"", lf[which(varNames == vars[i])], "\"/>", sep = "")
            cat(c("\n","\t","\t","<netcdf>"), sep = "", file = z)
            cat(c("\n","\t","\t","<aggregation dimName=\"time\" type=\"joinExisting\">"), sep = "", file = z)
            for (j in 1:length(subset)) {
                  cat(c("\n","\t","\t","\t",subset[j]), sep = "", file=z)
            }
            cat(c("\n","\t","\t","</aggregation>"), sep = "", file = z)
            cat(c("\n","\t","\t","</netcdf>"), sep = "", file = z)
      }
      cat(c("\n","\t","</aggregation>"), sep = "", file = z)
      cat(c("\n","</netcdf>"), sep = "", file = z)
      close(z)
      message(paste("[",Sys.time(),"]", sep=""))
      message(paste("NcML file \"",ncml.file, "\" created from ", length(lf), " files corresponding to ", length(vars), " variables", sep=""))
      message("Use 'dataInventory(NcML file)' to obtain a description of the dataset")
}
# End

# # Example 1: non-recursive dataset creation with filtering by years
# pat = paste(2001:2005, collapse="|")
# path <- "/home/juaco/.gvfs/SFTP for juaco on ui01.macc.unican.es/vols/seal/oceano/gmeteo/DATA/METEOLAB_PUBLIC/ObservationsData/WFDEI/daily/data/fwi"
# makeDataset(path, ncml.file= "../../exampleDatasets/fwi_WFDEI_daily.ncml", pattern=pat)
# # Example 2: recursive dataset creation with filtering by years
# pat = paste(2001:2005, collapse="|")
# path <- "/home/juaco/.gvfs/SFTP for juaco on ui01.macc.unican.es/vols/seal/oceano/gmeteo/DATA/METEOLAB_PUBLIC/ObservationsData/WFDEI/daily/data"
# makeDataset(source.dir=path, ncml.file="../../exampleDatasets/WFDEI_daily_2001_2005.ncml", pattern = pat, recursive = TRUE)
# # Example 3: recursive dataset creation with filtering by years and variables
# dataInventory("../../exampleDatasets/WFDEI_daily_2001_2005.ncml")
## [2013-12-08 15:07:39]
## Creating dataset from 497 files
## [2013-12-08 15:35:40]
## NcML file "../../exampleDatasets/WFDEI_daily_2001_2005.ncml" created from 497 files corresponding to 5 variables
## Use 'dataInventory(NcML file)' to obtain a description of the dataset