#' Finds vertical level from variable definition
#' 
#' @description Parses the variable name as passed by the 'var' argument and extracts the level and variable name
#' 
#' @param var Character string defining the (standandar) variable defined
#' @return A list with the variable name (string) and vertical level (double)
#' @details The level output of this function is passed to \code{\link{getVerticalLevelPars}}
#' @author J Bedia \email{joaquin.bedia@@gmail.com}
#' @keywords internal
#' @export

findVerticalLevel <- function(var) {
      if (grepl("@", var)) {
            aux <- unlist(strsplit(var, split = "@"))
            level <- as.double(aux[2])
            var <- aux[1]
      } else {
            level <- NULL
            var <- var
      }
      return(list("var" = var, "level" = level))
}
# End