makeDataset <- function(source.dir, ncml.file, pattern = NULL, recursive = FALSE) {
    file.create(ncml.file)
    if (is.null(pattern)) {
        lf <- normalizePath(list.files(source.dir, full.names = TRUE, pattern = "\\.nc$", recursive = recursive), winslash = "/")      
    } else {
        lf <- normalizePath(list.files(source.dir, full.names = TRUE, pattern = pattern, recursive = recursive), winslash = "/")      
        lf <- lf[grep("\\.nc$", lf)]
    }
    message("[",Sys.time(),"] Creating dataset from ", length(lf), " files")
    z <- file(ncml.file, "w")
    # http://www.unidata.ucar.edu/software/netcdf/ncml/v2.2/Aggregation.html
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
    if (identical(vars, varNames)) { # Unique files for each variable
        for (i in 1:length(vars)) {
            cat(c("\n","\t","<netcdf location=", "\"", lf[i]), "\"/>", sep = "", file = z)
        }
    } else { # Many files per variable
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
    }
    cat(c("\n","\t","</aggregation>"), sep = "", file = z)
    cat(c("\n","</netcdf>"), sep = "", file = z)
    close(z)
    message("[",Sys.time(),"] NcML file \"", ncml.file, "\" created from ", length(lf), " files corresponding to ", length(vars), " variables")
    message("Use 'dataInventory' to obtain a description of the dataset")
}
# End
# # dataset
# makeDataset(source.dir="/home/juaco/.gvfs/SFTP for juaco on ui01.macc.unican.es/vols/seal/oceano/gmeteo/DATA/ECA/EOBS/Grid_0.25deg_reg/v9/data/", ncml.file="ignore/prueba.ncml", recursive=TRUE, pattern="rr|tn|tx")
# di <- dataInventory("ignore//prueba.ncml")
# str(di)