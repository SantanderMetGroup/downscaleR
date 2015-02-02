#' @description Creates virtual datasets by modifying and combining other datasets via NcML.
#' @title Create a dataset from a collection of (netCDF) files
#' 
#' @import rJava

#' @details The NetCDF Markup Language (NcML) is an XML dialect that allows creating
#' CDM datasets (i.e.: any collection of scientific data which can be accessed
#' through the NetCDF-Java / CDM library). The NcML document refers to
#' another dataset called the referenced CDM dataset, generally composed on a
#' number of netCDF files (but also grib, hdf or many other binary file formats)
#' containing the geo-referenced data. This function creates a NcML file from
#' multiple CDM files that are conveniently combined ("aggregated") along the
#' selected dimension.
#' The use of NcML is not only intended for CDM file combination, and this function considers just a particular case.
#' NcML is a powerful, yet relatively simple way of dealing with large, complex datasets in a straightforward manner.
#' Among other capabilitities, NcML files can be generated in order to add/delete
#' metadata and variables to be renamed, added, deleted and restructured.
#' This function is intended for a simple operation of aggregation of collections of netCDF files,
#' as it is the most common case of gridded climate datasets.
#'
#' @param source.dir Parent directory containing the files to be aggregated
#' @param ncml.file Full path of the output NcML file
#' @param file.ext Character string indicating the extension of the CDM datasets
#'  to be aggregated. Default to \code{nc} (netCDF).
#' @param aggr.dim Character string indicating the dimension along which the 
#' files will be concatenated. Default to \code{"time"}.
#' @param pattern An optional regular expression. Only file names which match the
#' regular expression will be considered in the aggregation 
#' (see \link[base]{regexp}). This argument can be useful in order to save time
#' when only a particular subset of variables from the whole collection is 
#' needed. Default to \code{NULL}, meaning that all files in the search path
#' are included (See next argument).
#' @param recursive Logical. Should the listing of files to be aggregated recurse
#' into directories?. Default to \code{FALSE}. This is useful for instance when 
#' each variable is stored in a sepparate subdirectory.
#' @param verbose Logical. Should additional information of the NcML file creation
#' steps be printed on screen?. Default to \code{TRUE}.
#' @return Creates a NcML file at the specified location.
#' @export
#' @note The current implementation of the function only considers datasets in which each file stores
#' one single variable. For other dataset configurations, please refer to the NcML tutorial.
#' @references NcML Tutorial \url{http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/ncml/Tutorial.html}. Accessed 28 Feb 2015.).
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @family loading

makeAggregatedDataset <- function(source.dir, ncml.file, file.ext = "nc", aggr.dim = "time",
                        pattern = NULL, recursive = FALSE, verbose = TRUE) {
      file.ext <- match.arg(file.ext, c("nc", "hdf", "grib","gini"))
      suffix <- paste("\\.", file.ext, "$", sep = "")
      if (is.null(pattern)) {
            pattern <- suffix
      }
      lf <- normalizePath(list.files(source.dir, full.names = TRUE, pattern = pattern, recursive = recursive), winslash = "/")
      lf <- grep(suffix, lf, value = TRUE)
      message("[",Sys.time(),"] Creating dataset from ", length(lf), " files")
      z <- file(ncml.file, "w")
      # http://www.unidata.ucar.edu/software/netcdf/ncml/v2.2/Aggregation.html
      cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file = z)
      cat(c("\n","<netcdf xmlns=\"http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2\">"), sep = "", file = z)
      cat(c("\n","\t","<aggregation type=\"union\">"), sep = "", file = z)
      varNames <- rep(NA, length(lf))
      for (i in 1:length(lf)) {
            gds <- J("ucar.nc2.dt.grid.GridDataset")$open(.jnew("java/lang/String", lf[i]))
            varNames[i] <- unlist(strsplit(gsub("\\[|]|\\s", "", gds$getGrids()$toString()), ","))[1] # case when you have something like "varName, lon, lat" (e.g. some ENSEMBLES files)
      }
      vars <- unique(varNames)
      #¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
      # Next code chunk takes one unique file for each variable and checks
      # if the aggregation dimension exists for that variable.
      # This makes sense when aggregation e.g. by "level", and we want to include
      # surface variables too in the dataset
      #¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
      hasAggrDim <- rep(TRUE, length(vars))
      ind <- match(vars, varNames)
      for (i in 1:length(vars)) {
            if (verbose) {
                  message("[",Sys.time(),"] Scanning file ", i, " out of ", length(vars))
            }
            gds <- J("ucar.nc2.dt.grid.GridDataset")$open(lf[ind[i]])
            grid <- gds$findGridByShortName(vars[i])
            dims <- unlist(strsplit(gsub("\\[|]|;|\\s","", grid$getDimensions()$toString()), split = ","))
            hasAggrDim[i] <- any(grepl(aggr.dim, dims))
            gds$close()
      }
      #¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
      # Those variables without the aggregating dimension are directly aggregated
      # with "union" at the upper level
      #¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
      ind <- which(!hasAggrDim)
      if (length(ind) > 0) {
            for (i in 1:length(ind)) {
                  varfile <- grep(vars[ind[i]], lf, value = TRUE)
                  for (j in 1:length(varfile)) {
                        cat(c("\n","\t","<netcdf location=", "\"", varfile[j]), "\"/>", sep = "", file = z)
                  }
            }
      }
      #¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
      # Case when there is a single file per variable. No need to concatenate
      # different files of the same variable. The variables are directly aggregated
      # with "union" at the upper level
      #¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
      if (length(lf) == length(vars)) {
            for (i in 1:length(lf)) {
                  gds <- J("ucar.nc2.dt.grid.GridDataset")$open(lf[i])
                  grid <- gds$findGridByShortName(vars[i])
                  dims <- unlist(strsplit(gsub("\\[|]|;|\\s","", grid$getDimensions()$toString()), split = ","))
                  aggrDimIndex <- as.integer(grep(aggr.dim, dims) - 1)
                  ncoords <- grid$getDimension(aggrDimIndex)$getLength()
                  cat(c("\n", "\t", "<netcdf location=", "\"", lf[i], "\" ncoords=\"", ncoords), "\"/>", sep = "", file = z)
            }
      } else {
      #¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
      # Many files per varible are concatenated by "joinExisting" at the next
      # aggregation level
      #¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
            ind <- which(hasAggrDim)
            if (length(ind) < 1) {
                  stop("Invalid dimension for aggregation")
            }
            if (verbose) {
                  message("[",Sys.time(),"] Defining aggregating dimension length\nThis process may be slow but will significantly speed-up data retrieval...")
            }
            for (i in 1:length(ind)) {
                  varfile <- grep(vars[ind[i]], lf, value = TRUE)
                  cat(c("\n","\t","\t","<netcdf>"), sep = "", file = z)
                  cat(c("\n","\t","\t","<aggregation dimName=\"", aggr.dim, "\" ", "type=\"joinExisting\">"), sep = "", file = z)
                  cat(c("\n","\t","\t","<variableAgg name=\"", vars[ind[i]], "\"/>"), sep = "", file = z)
                  for (j in 1:length(varfile)) {
                        gds <- J("ucar.nc2.dt.grid.GridDataset")$open(varfile[j])
                        grid <- gds$findGridByShortName(vars[ind[i]])
                        dims <- unlist(strsplit(gsub("\\[|]|;|\\s","", grid$getDimensions()$toString()), split = ","))
                        aggrDimIndex <- as.integer(grep(aggr.dim, dims) - 1)
                        ncoords <- grid$getDimension(aggrDimIndex)$getLength()
                        gds$close()
                        cat(c("\n","\t","\t","\t","<netcdf location=", "\"", varfile[j]), "\" ncoords=\"", ncoords, "\"/>", sep = "", file = z)
                  }
                  cat(c("\n","\t","\t","</aggregation>"), sep = "", file = z)
                  cat(c("\n","\t","\t","</netcdf>"), sep = "", file = z)
            }
            if (verbose) {
                  message("[",Sys.time(),"] Dimension length defined")
            }
      }
      cat(c("\n","\t","</aggregation>"), sep = "", file = z)
      cat(c("\n","</netcdf>"), sep = "", file = z)
      close(z)
      message("[",Sys.time(),"] NcML file \"", ncml.file, "\" created from ", length(lf), " files corresponding to ", length(vars), " variables")
      message("Use 'dataInventory' to obtain a description of the dataset")
}
# End


