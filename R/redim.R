#' @title Complete missing dimensions of Grid objects
#' @description Inverse of drop to complete all dimensions of the Data array
#' @param grid A grid 
#' @return The same object with all the dimensions (i.e. member, time, station)
#' @keywords internal
#' @importFrom abind abind
#' @author M. Iturbide

redim <- function(obj, drop = FALSE) {
      if (drop == FALSE){
            if("station" %in% attr(obj$Data, "dimensions")){
                  ind <- which(attr(obj$Data, "dimensions")=="station")
                  dimNames <- c(attr(obj$Data, "dimensions")[-ind], "lat", "lon")
                  obj$Data <- unname(abind(obj$Data, NULL, along = ind+1))
                  attr(obj$Data, "dimensions") <- dimNames
            }else{
            # Add fake 'member' dimension to single-station datasets
                  if (!("member" %in% attr(obj$Data, "dimensions"))) {
                        dimNames <- c("member", attr(obj$Data, "dimensions"))
                        obj$Data <- unname(abind(obj$Data, NULL, along = 0))
                        attr(obj$Data, "dimensions") <- dimNames
                  }
            # Add fake runtime dimension to deterministic/obs
                  if (!("runtime" %in% attr(obj$Data, "dimensions"))) {
                        dimNames <- c("runtime", attr(obj$Data, "dimensions"))
                        obj$Data <- unname(abind(obj$Data, NULL, along = -1))    
                        attr(obj$Data, "dimensions") <- dimNames
                  }
            }
      }else{
            if(1 %in% dim(obj$Data)){
                  dimNames <- attr(obj$Data, "dimensions")[-which(dim(obj$Data) == 1)]
                  obj$Data <- drop(obj$Data)
                  attr(obj$Data, "dimensions") <- dimNames
                  if("lat" %in% dimNames & "lon" %in% dimNames == FALSE)
                        attr(obj$Data, "dimensions")[attr(obj$Data, "dimensions") == "lat"] <- "station"
            }else{
                  obj <- obj
            }
            
      }
      return(obj) 
            
}
         
#End

