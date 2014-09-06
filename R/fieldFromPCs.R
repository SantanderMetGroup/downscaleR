#' @title Reconstruct a Field from EOFs and principal components
#' 
#' @description Reconstructs a field of a climatic variable from the outputs of a
#' principal components analysis
#' 
#' @param prinCompObj A EOF analysis object as returned by \code{\link{prinComp}}
#' @param var Character string indicating the variable to be re-constructed. In case of
#' PCA analyses performed on a single variable (either multimember or not), this can be omitted.
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
#' @seealso \code{\link{prinComp}}, \code{\link{plotEOF}}
#' 
#' @examples \dontrun{
#' # First a multifield containing a set of variables is loaded (e.g. data for spring spanning the 30-year period 1981--2010):
#' ncep <- file.path(find.package("downscaleR"), "datasets/reanalysis/Iberia_NCEP/Iberia_NCEP.ncml")
#' multifield <- loadMultiField(ncep, vars = c("hus@@85000", "ta@@85000", "psl"), season = c(3:5), years = 1981:2010)
#' # In this example, we retain the first 10 PCs
#' pca <- prinComp(multifield, n.eofs = 10)
#' # We now recover the sea-level pressure field from the PCs:
#' names(pca)
#' psl2 <- fieldFromPCs(pca, "psl")
#' str(psl2)
#' # The attributes of psl2 indicate that this is a reconstructed field from 10 PCs, explaining 99\% of the variance:
#' attributes(psl2)
#' multifield$Variable$varName
#' # psl is the 3rd one
#' # The mean fields of both the original and the reconstructed fields is computed:
#' psl.orig <- multifield$Data[3,,,]
#' psl.reconstructed <- psl2$Data
#' z <- apply(psl.orig, c(3,2), mean)
#' z1 <- apply(psl.reconstructed, c(3,2), mean)
#' # These are the spatial coordinates
#' x <- psl2$xyCoords$x
#' y <- psl2$xyCoords$y
#' require(fields)
#' par(mfrow = c(2,2))
#' image.plot(x,y,z, asp = 1, horizontal = TRUE)
#' world(add = TRUE)
#' title("Original SLP field")
#' image.plot(x,y,z1, asp = 1, horizontal = TRUE)
#' world(add = TRUE)
#' title("Reconstructed SLP field")
#' mtext("(Using the first 10 PCs)")
#' image.plot(x,y,z1-z, asp = 1)
#' world(add = TRUE)
#' title("Difference (bias)")
#' par(mfrow = c(1,1))
#'  
#' # An example of multimember reconstruction from a multimember PC analysis:
#' data(tasmax_forecast)
#' # Note that multimember pca analysis takes some time, depending on the computer
#' pca2 <- prinComp(tasmax_forecast, n.eofs = 10)
#' tasmax_recovered <- fieldFromPCs(pca2)
#' plotMeanField(tasmax_recovered, multi.member = TRUE)
#' # Also note that now the length of the "nPCs" and "explained_variance" matches the number of members
#' attributes(tasmax_recovered)
#' }

fieldFromPCs <- function(prinCompObj, var) {
      varNames <- attributes(prinCompObj)$names 
      if (length(varNames) == 1L) {
            var <- varNames
            var.ind <- 1
      } else {
            if (var == "COMBINED") {
                  stop("The combined field can't be reconstructed (only individual variables)")
            }
            var.ind <- match(var, varNames)
            if(is.na(var.ind)) {
                  stop("Variable not found.\nCheck the 'variables' attribute of the input")
            }
      }
      scaling <- attributes(prinCompObj)$"scaled:method"
      x <- attributes(prinCompObj)$xCoords
      y <- attributes(prinCompObj)$yCoords
      pco <- prinCompObj[[var.ind]]
      prinCompObj <- NULL
      scale <- lapply(1:length(pco), function(x) attributes(pco[[x]])$"scaled:scale")
      center <- lapply(1:length(pco), function(x) attributes(pco[[x]])$"scaled:center")
      exv <- lapply(1:length(pco), function(x) attributes(pco[[x]])$"explained_variance")
      level <- attributes(pco)$level  
      Members <- names(pco)
      aux.list  <- lapply(1:length(pco), function(n) {
            Xhat <- (pco[[n]]$PCs %*% t(pco[[n]]$EOFs)) * scale[[n]] + center[[n]]
            out <- mat2Dto3Darray(Xhat, x, y)
            return(out)
      })
      pco <- NULL
      dimNames <- attr(aux.list[[1]], "dimensions")
      Data <- unname(do.call("abind", c(aux.list, along = -1)))
      if (identical(Data, drop(Data))) {
            dimNames <- append(dimNames, "member", after = 0)
            mm <- TRUE
      } else {
            Data <- drop(Data)
            mm <- FALSE
      }
      attr(Data, "dimensions") <- dimNames
      out <- list("Variable" = list("varName" = var, "level" = level), "Data" = Data, "xyCoords" = list("x" = x, "y" = y), "Members" = Members)
      if (!mm) {
            out <- out[-length(out)]      
      }
      attr(out, "nPcs") <- unlist(lapply(1:length(exv), function(x) length(exv[[x]])))
      attr(out, "explained_variance") <- unlist(lapply(1:length(exv), function(x) round(tail(exv[[x]], 1), 2)))
      return(out)
}
# End




