#' Finds vertical level from variable definition
#' 
#' @param var Character string defining the (standandar) variable defined
#' @return A numeric vector of length one with the vertical level defined
#' @details The output of this function is passed to \code{\link{getVerticalLevelPars}}
#' @author J Bedia \email{joaquin.bedia@@gmail.com}

findVerticalLevel <- function(var) {
      if (grepl("@", var)) {
            level <- tail(unlist(strsplit(var, split = "@")), 1)
      } else {
            level <- NULL
      }
      return(level)
}
# End