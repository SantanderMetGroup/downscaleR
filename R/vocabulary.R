#' @title Standard definition of climate variables
#' 
#' @description Definition of a standard naming convention (and units) for climate variables,
#' with the following elements:
#' 
#' \itemize{
#'   \item identifier. R standard name of the climate variable
#'   \item standard_name. description of the variable
#'   \item units. Standard variable units
#' }
#' 
#' @family homogenization
#' @details The dictionary provides the reference names for variable homogenization. All the required
#' variable transformations to convert the original variable (as returned by \code{\link{dataInventory}}) and
#' the standard variables as denoted in the vocabulary's \code{identifier} are specified by simple parameters
#' in the \file{dictionary} (.dic) file. 
#' 
#' @template templateDicDetails
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with \emph{n} rows (the number of variables may change with the package 
#' version as new variables are included) and 3 variables.
#' @name vocabulary
#' @examples
#' data(vocabulary)
#' print(vocabulary)
NULL