#' @title Reconstruct a Field from EOFs and principal components
#' 
#' @description Reconstructs a field of a climatic variable from the outputs of a
#' principal components analysis
#' 
#' @param prinCompObj A EOF analysis object as returned by \code{\link{prinComp}}
#' @param var Character string indicating the variable to be re-constructed
#' 
#' @return A list similar to the \code{\link{loadGridData}} output, but simplified. See details.
#' 
#' @details The output of this function returns the minimum required information to use the
#'  \code{\link{plotMeanField}} method, and is intended for comparison and visual analysis
#'  of the differences between the original fields and the reconstructed ones, for instance
#'  in determining an optimal number of PCs etc... Hence, the temporal information (i.e.,
#'   the \code{Dates} object) is lost, and should be retrieved from the original field/multifield
#'   object used to compute the PC/EOF analysis.
#' 
#' @importFrom abind abind
#' 
#' @export
#' 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' 
#' @family loading.grid
#' @family multifield
#'  


# prinCompObj <- prinComp(gridData.mf, v.exp = .95)
# str(prinCompObj)
# var = "q850mb"

fieldFromPCs <- function(prinCompObj, var) {
      varNames <- attributes(prinCompObj)$variables #[1:n.vars]
      if (var == "COMBINED") {
            stop("The combined field cannot be reconstructed (only individual variables)")
      }
      var.ind <- match(var, varNames)
      if(is.na(var.ind)) {
            stop("Variable not found.\nCheck the 'variables' attribute of the input")
      }
      scaling <- attributes(prinCompObj)$"scaled:method"
      x <- attributes(prinCompObj)$xCoords
      y <- attributes(prinCompObj)$yCoords
      pco <- prinCompObj[[var.ind]]
      prinCompObj <- NULL
      scale <- attributes(pco)$"scaled:scale"
      center <- attributes(pco)$"scaled:center"
      exv <- attributes(pco)$explained_variance
      Xhat <- (pco$PCs %*% t(pco$EOFs)) * scale + center
      Data <- mat2Dto3Darray(Xhat, x, y)
      out <- list("Variable" = list("varName" = var), "Data" = Data, "xyCoords" = list("x" = x, "y" = y))
      attr(out, "nPCs") <- length(exv)
      attr(out, "explained_variance") <- round(tail(exv, 1),2)
      return(out)
}

# a <- fieldFromPCs(prinComp(gridData.mf, v.exp = .99), "q850mb")
# str(a)
# plotMeanField(gridData.mf)
# gridData <- a

# plotMeanField(gridData.mf)

