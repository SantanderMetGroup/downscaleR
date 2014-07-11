#' Inventory of a gridded dataset
#' 
#' Returns a list with summary information about the variables stored in a gridded dataset.
#' Sub-routine of \code{dataInventory}
#' 
#' @param dataset A full path to the file describing the dataset (NcML)
#' @return A (named) list whose length is determined by the number of variables stored in the dataset,
#' its names corresponding to the short names of the variables.
#' For each variable, information on the variable long name, data type, units and
#' characteristics of its dimensions is provided.
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @keywords internal
 

dataInventory.NetCDF <- function(dataset) {
    gds <- J("ucar.nc2.dt.grid.GridDataset")$open(dataset)
    varNames <- unlist(strsplit(gsub("\\[|]|\\s", "", gds$getGrids()$toString()), ","))
    rm.ind <- grep("^lon|^lat", varNames)  
    if (length(rm.ind) > 0) {
        varNames <- varNames[-rm.ind]
    }
    if (length(varNames) == 0) {
        stop("No variables found")
    } else {
        var.list <- list()
        for (i in 1:length(varNames)) {
            message("[", Sys.time(), "] Retrieving info for \'", varNames[i], "\' (", length(varNames) - i, " vars remaining)")
            description <- gds$getDataVariable(varNames[i])$getDescription()
            varName <- gds$getDataVariable(varNames[i])$getShortName()
            dataType <- gds$getDataVariable(varNames[i])$getDataType()$toString()
            units <- gds$getDataVariable(varNames[i])$getUnitsString()
            grid <- gds$findGridByShortName(varName)
            dim.list <- scanVarDimensions(grid)
            var.list[[i]] <- list("Description" = description, "DataType" = dataType, "Units" = units, "Dimensions" = dim.list)
        }
        names(var.list) <- varNames
    }
    gds$close()
    return(var.list)
}
# End