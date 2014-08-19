#' @title Plot an arbitrary number of EOFs
#' 
#' @description Plots an arbitrary number of EOFs. Useful to have a quick overview of the main spatial modes
#'  of a field.
#' @param prinCompObj A PCA object as returned by \code{\link{prinComp}}
#' @param var Character string indicating the variable whose EOFs are to be displayed. If the PCA analysis has
#' been applied to 1 single field, this argument can be omitted.
#' @param n.eofs Number of EOFs to be displayed. Default to NULL, indicating that all computed EOFS
#' will be represented
#' 
#' @return A plot with as many panels as EOFs requested, in the original units of the variable
#' @export
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' 
#' @family multifield
#' @family loading.grid
#' 
#' @examples \dontrun{ 
#' # Winter temperature at 850 mb isobaric surface pressure level is loaded (period 1981--2010):
#' ncep <- file.path(find.package("downscaleR"), "datasets/reanalysis/Iberia_NCEP/Iberia_NCEP.ncml")
#' field <- loadGridData(ncep, var = "ta@@85000", season = c(12,1,2), years = 1981:2010)
#' # PCA analysis, retaining the PCs that explain 95\% of the total variance:
#' pca <- prinComp(field, v.exp = .95)
#' # Plot of all EOFs
#' plotEOF(pca)
#' # Plot the first 4 EOFs:
#' plotEOF(pca, n.eofs = 4)
#' 
#' # Example with PCA analysis of a multifield (multiple variables)
#' load geopotential heigth at 500 mb, temperature at 1000 mb and sea-level pressure
#' multifield <- loadMultiField(ncep, vars = c("z@@50000", "ta@@100000", "psl"), season = c(12,1,2), years = 1981:2010)
#' # PCA analysis, retaining the first 9 PCs of each variable:
#' pca2 <- prinComp(multifield, n.eofs = 9)
#' names(pca2)
#' # EOFs for geopotential
#' plotEOF(pca2, "z")
#' plotEOF(pca2, "ta", n.eofs = 4)
#' plotEOF(pca2, "psl", n.eofs = 2)
#' }


plotEOF <- function(prinCompObj, var = NULL, n.eofs = NULL) {
      varNames <- attributes(prinCompObj)$names
      if (length(varNames) == 1) {
            ind.var <- 1
      } else {
            if (is.null(var)) {
                  stop("The argument 'var' is missing, with no default value")
            }
            ind.var <- match(var, varNames)
            if (is.na(ind.var)) {
                  stop("Variable given in 'var' not found")
            }
            if (var == "COMBINED") {
                  stop("It is not possible to display the combined EOF")
            }
      }
      pco <- prinCompObj[[ind.var]]
      x <- attributes(prinCompObj)$xCoords
      y <- attributes(prinCompObj)$yCoords
      if (is.null(n.eofs)) {
            n.eofs <- ncol(pco$EOFs)
      } else if (n.eofs > ncol(pco$EOFs)) {
            warning("The requested number of EOFs (",n.eofs,") is larger than actual number of EOFs available (",ncol(pco$EOFs),")\nOnly the available EOFs were displayed")
            n.eofs <- ncol(pco$EOFs)
      }
      mu <- attributes(pco)$"scaled:center"
      sigma <- attributes(pco)$"scaled:scale"
      out <- list("Data" = mat2Dto3Darray(t(pco$EOFs[ ,1:n.eofs] * sigma + mu), x, y), "xyCoords" = list(x = x, y = y))
      multiPlot(out, split.dim.name = "time", title = paste("EOF",1:n.eofs))
}
# End
