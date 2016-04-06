    ## grid2PCs.R Projection of a grid onto an EOF

    ## Copyright (C) 2016 Santander Meteorology Group (http://www.meteo.unican.es)

    ## This program is free software: you can redistribute it and/or modify
    ## it under the terms of the GNU General Public License as published by
    ## the Free Software Foundation, either version 3 of the License, or
    ## (at your option) any later version.

    ## This program is distributed in the hope that it will be useful,
    ## but WITHOUT ANY WARRANTY; without even the implied warranty of
    ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    ## GNU General Public License for more details.

    ## You should have received a copy of the GNU General Public License
    ## along with this program.  If not, see <http://www.gnu.org/licenses/>. 


#' @title Projection of a grid onto an EOF
#' @description Projection of a grid onto a user-defined EOF
#' @param prinCompObj An object created by \code{\link{prinComp}}
#' @param grid Input grid to project. Multigrids are not allowed.
#' @param n.pcs Number of principal components to be retained. Default to the number of EOFs contained in the input \code{prinCompObj}.
#' @return A list, each component corresponding to a member (a list of length 1 if the projected grid is not multimember),
#' containing a matrix of PCs, being the number of columns defined by \code{n.pcs}, and the number of rows by the
#' length of the time dimension.
#' @details The function is intended to project a grid onto user specified EOFs, contained in \code{prinCompObj}.
#' Note that the grid should have the same spatial extent and resolution than the original EOF, so input grids from
#' different models should be first interpolated. Also, to ensure consistency, both objects should contain the same 
#' variable name.  
#' 
#' Note that the function is not currently implemented to deal with decadal predictions
#'  (i.e., 'runtime' dimension is not handled)
#' @author J Bedia
#' @export
#' @seealso \code{\link{prinComp}} for EOF analysis.


grid2PCs <- function(prinCompObj, grid, n.pcs = NULL) {
      gridName <- grid$Variable$varName
      if (!(gridName %in% names(prinCompObj))) {
            stop("Input PCA object and grid name do not match", call. = FALSE)
      }
      EOF <- prinCompObj[[gridName]][[1]][["EOFs"]]
      mu <- attr(prinCompObj[[gridName]][[1]], "scaled:center")
      sigma <- attr(prinCompObj[[gridName]][[1]], "scaled:scale")
      prinCompObj <- NULL
      grid <- redim(grid)
      dimNames <- attr(grid$Data, "dimensions")
      lat.ind <- grep("lat", dimNames)
      lon.ind <- grep("lon", dimNames)
      if (prod(dim(grid$Data)[c(lat.ind,lon.ind)]) != nrow(EOF)) {
            stop("Incompatible array dimensions. Input grid and EOF must be in the same grid")
      }
      if (is.null(n.pcs)) {
            n.pcs <- ncol(EOF)
      } else if (n.pcs > ncol(EOF)) {
            message("Number of PCs requested exceeds EOF dimensions. ", ncol(EOF), " PCs will be returned.")
            n.pcs <- ncol(EOF)
      }
      mem.ind <- grep("member", dimNames)
      n.mem <- dim(grid$Data)[mem.ind]
      rt.ind <- grep("runtime", dimNames)
      n.rt <- dim(grid$Data)[rt.ind]
      for (i in 1:n.rt) {# Runtimes are currently ignored, and assumed to be only 1
            PC.list <- lapply(1:n.mem, function(j) {
                  X <- grid$Data[i,j,,,]
                  attr(X, "dimensions") <- c("time","lat","lon")
                  X <- array3Dto2Dmat(X)
                  X <- (X - mu) / sigma
                  PCs <- X %*% EOF 
                  PCs[,1:n.pcs,drop = FALSE]
            })
      }
      return(PC.list)
}
# End

