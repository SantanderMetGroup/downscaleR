#' @title Complete missing dimensions of Grid or Station objects
#' @description Complete all dimensions of the Data array
#' @param obj A grid or station data
#' @param runtime logical. Add runtime dimension (default = FALSE)
#' @param drop logical. Drop dimensions of length = 1 (default = FALSE)
#' @return The same object with all the dimensions (i.e. member, time, station)
#' @keywords internal
#' @importFrom abind abind
#' @author M. Iturbide

redim <- function(obj, runtime = FALSE, drop = FALSE) {
      if (drop == FALSE) {
            if ("station" %in% attr(obj$Data, "dimensions")) {
                  ind <- which(attr(obj$Data, "dimensions") == "station")
                  dimNames <- c(attr(obj$Data, "dimensions")[-ind], "lat", "lon")
                  obj$Data <- unname(abind(obj$Data, NULL, along = ind + 1))
                  attr(obj$Data, "dimensions") <- dimNames
            } else {
                  # Add fake 'time' dimension  -------
                  if (!("time" %in% attr(obj$Data, "dimensions"))) {
                        dimNames <- c("time", attr(obj$Data, "dimensions"))
                        obj$Data <- unname(abind(obj$Data, NULL, along = 0))
                        attr(obj$Data, "dimensions") <- dimNames
                  }
                  # Add fake 'member' dimension  -------
                  if (!("member" %in% attr(obj$Data, "dimensions"))) {
                        dimNames <- c("member", attr(obj$Data, "dimensions"))
                        obj$Data <- unname(abind(obj$Data, NULL, along = 0))
                        attr(obj$Data, "dimensions") <- dimNames
                  }
                  # Add fake runtime dimension to deterministic/obs -----------
                  if (!("runtime" %in% attr(obj$Data, "dimensions")) & runtime == TRUE) {
                        dimNames <- c("runtime", attr(obj$Data, "dimensions"))
                        obj$Data <- unname(abind(obj$Data, NULL, along = -1))    
                        attr(obj$Data, "dimensions") <- dimNames
                  }
            }
            dimNames <- c("runtime", "member","time","lat","lon")
            dimNames.aux <- attr(obj$Data, "dimensions")
            perm <- na.omit(match(dimNames, dimNames.aux))
            obj$Data <- aperm(obj$Data, perm)
            obj[["Data"]] <- unname(obj$Data)
            attr(obj[["Data"]], "dimensions") <- dimNames.aux[perm]
      } else {
            if (1 %in% dim(obj$Data)) {
                  dimNames <- attr(obj$Data, "dimensions")[-which(dim(obj$Data) == 1)]
                  obj$Data <- drop(obj$Data)
                  attr(obj$Data, "dimensions") <- dimNames
                  if ("lat" %in% dimNames & "lon" %in% dimNames == FALSE) {
                        attr(obj$Data, "dimensions")[attr(obj$Data, "dimensions") == "lat"] <- "station"
                  }
            }
      }
      return(obj) 
}
#End

