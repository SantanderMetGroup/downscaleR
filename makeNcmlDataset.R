makeNcmlDataset <- function(source.dir, ncml.file) {
      source.dir = source.dir
      ncml.file = ncml.file
      file.create(ncml.file)
      lf <- normalizePath(list.files(source.dir, full.names=TRUE, pattern = "\\.nc$"), winslash = "/")
      z <- file(ncml.file, "w")
      # [http://www.unidata.ucar.edu/software/netcdf/ncml/v2.2/Aggregation.html]
      cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", file = z)
      cat(c("\n","<netcdf xmlns=\"http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2\">"), sep = "", file = z)
      cat(c("\n","\t","<aggregation type=\"union\">"), sep = "", file = z)
      varNames <- rep(NA, length(lf))
      for (i in 1:length(lf)) {
            jString = .jnew("java/lang/String", lf[i])
            gds = J("ucar.nc2.dt.grid.GridDataset")$open(jString)
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
