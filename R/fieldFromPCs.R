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
#' @examples \dontrun{
#' # First a multifield containing a set of variables is loaded (e.g. data for spring spanning the 30-year period 1981--2010):
#' ncep <- file.path(find.package("downscaleR"), "datasets/reanalysis/Iberia_NCEP/Iberia_NCEP.ncml")
#' multifield <- loadMultiField(ncep, vars = c("hus@85000", "ta@85000", "psl"), season = c(3:5), years = 1981:2010)
#' # In this example, we retain the first 10 PCs
#' pca <- prinComp(multifield, n.eofs = 10)
#' # We now recover the sea-level pressure field from the PCs:
#' names(pca)
#' psl2 <- fieldFromPCs(pca, "psl")
#' str(psl2)
#' # The attributes of psl2 indicate that this is a reconstructed field from 10 PCs, explaining 99\% of the variance:
#' attributes(psl2)
#'multifield$Variable$varName
#'# psl is the 3rd one
#'# The mean fields of both the original and the reconstructed fields is computed:
#'psl.orig <- multifield$Data[3,,,]
#'psl.reconstructed <- psl2$Data
#'z <- apply(psl.orig, c(3,2), mean)
#'z1 <- apply(psl.reconstructed, c(3,2), mean)
#'# These are the spatial coordinates
#'x <- psl2$xyCoords$x
#'y <- psl2$xyCoords$y
#'par(mfrow = c(2,2))
#'image.plot(x,y,z, asp = 1, horizontal = TRUE)
#'world(add = TRUE)
#'title("Original SLP field")
#'image.plot(x,y,z1, asp = 1, horizontal = TRUE)
#'world(add = TRUE)
#'title("Reconstructed SLP field")
#'mtext("(Using the first 10 PCs)")
#'image.plot(x,y,z1-z, asp = 1)
#'world(add = TRUE)
#'title("Difference (bias)")
#'par(mfrow = c(1,1))
#'}
#'

fieldFromPCs <- function(prinCompObj, var) {
      varNames <- attributes(prinCompObj)$names #[1:n.vars]
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
# End

