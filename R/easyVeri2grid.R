#' @title easyVerification matrix to climatological grid conversion
#' @description Convert a xyz-type score matrix as returned by \code{veriApply} to a climatological grid
#' @param easyVeri.mat A matrix containing the verification score,
#'  as returned by \code{\link[easyVerification]{veriApply}}
#' @param obs.grid The grid containing the verifying observations used in the call to \code{veryApply}
#'  producing the score matrix.
#' @param verifun Optional. Character string indicating the value of the \code{verifun} value. Just for 
#' a better traceability and metadata completeness. 
#' @return A climatological grid.
#' @seealso \code{\link{climatology}}, \code{\link{plotClimatology}}.
#' 
#' @export
#' 
#' @author J. Bedia

easyVeri2grid <- function(easyVeri.mat, obs.grid, verifun = NULL) {
      x <- obs.grid$xyCoords$x     
      y <- obs.grid$xyCoords$y     
      if (length(x) != ncol(easyVeri.mat) | length(y) != nrow(easyVeri.mat)) {
            stop("XY coordinates and matrix dimensions do not match")      
      }
      obs.grid$Data <- easyVeri.mat
      attr(obs.grid$Data, "dimensions") <- c("lat", "lon")
      obs.grid <- redim(obs.grid, member = FALSE)
      # Fake climatology:fun attribute
      clim.att <- ifelse(is.null(verifun), "easiVeri", verifun)
      attr(obs.grid[["Data"]], "climatology:fun") <- clim.att
      return(obs.grid)
}

