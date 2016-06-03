#' @title Complete missing dimensions of Grid or Station objects
#' @description Complete all dimensions of the Data array
#' @param obj A grid or station data
#' @param member logical. Add member dimension (default = TRUE)
#' @param runtime logical. Add runtime dimension (default = FALSE)
#' @param drop logical. Drop dimensions of length = 1 (default = FALSE)
#' @return The same object with all the dimensions (i.e. member, time, station)
#' @details The function down not handle multigrids (i.e. \code{"var"} dimension).
#'  Thus, multigrids need to be subsetted along \code{"var"} prior to redimensioning.
#' @keywords internal
#' @importFrom abind abind
#' @importFrom stats na.omit
#' @author M. Iturbide, J. Bedia

redim <- function(obj,
                  member = TRUE,
                  runtime = FALSE,
                  drop = FALSE) {
      stopifnot(is.logical(member) | is.logical(runtime) | is.logical(drop))
      dimNames <- getDim(obj)
      if (!drop) {
            if ("station" %in% dimNames) {
                  # Add fake 'coordinates dimension' dimension  -------
                  ind <- match("station", dimNames)
                  dimNames <- c(dimNames[-ind], "lat", "lon")
                  obj$Data <- unname(abind(obj$Data, NULL, along = ind + 1))
                  attr(obj$Data, "dimensions") <- dimNames
            }
            # Add fake 'time' dimension  -------
            if (!("time" %in% dimNames)) {
                  dimNames <- c("time", dimNames)
                  obj$Data <- unname(abind(obj$Data, NULL, along = 0))
                  attr(obj$Data, "dimensions") <- dimNames
            }
            # Add fake 'member' dimension  -------
            if (!("member" %in% dimNames) & isTRUE(member)) {
                  dimNames <- c("member", dimNames)
                  obj$Data <- unname(abind(obj$Data, NULL, along = 0))
                  attr(obj$Data, "dimensions") <- dimNames
            }
            # Add fake runtime dimension to deterministic/obs -----------
            if (!("runtime" %in% dimNames) & isTRUE(runtime)) {
                  dimNames <- c("runtime", dimNames)
                  obj$Data <- unname(abind(obj$Data, NULL, along = -1))    
                  attr(obj$Data, "dimensions") <- dimNames
            }
            dimNames <- c("runtime", "member", "time", "lat", "lon")
            dimNames.aux <- attr(obj$Data, "dimensions")
            perm <- na.omit(match(dimNames, dimNames.aux))
            obj$Data <- aperm(obj$Data, perm)
            obj[["Data"]] <- unname(obj$Data)
            attr(obj[["Data"]], "dimensions") <- dimNames.aux[perm]
      } else {
            shp <- dim(obj[["Data"]])
            if (1 %in% shp) {
                  dimNames <- dimNames[-match(1, shp)]
                  obj$Data <- drop(obj$Data)
                  attr(obj$Data, "dimensions") <- dimNames
                  if ("lat" %in% dimNames & !("lon" %in% dimNames)) {
                        attr(obj$Data, "dimensions")[attr(obj$Data, "dimensions") == "lat"] <- "station"
                  }
            }
      }
      return(obj) 
}
#End

