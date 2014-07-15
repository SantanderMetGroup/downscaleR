#' @return A list with the following elements providing the necessary information 
#' for data representation and analysis:
#' \item{\code{Variable }}{A list with three elements:}
#'      \itemize{ 
#'            \item \code{varName} A character string indicating which is the variable returned. 
#'            Same as value provided for argument \code{var}
#'            \item \code{isStandard} Logical value indicating whether the variable returned 
#'            is standard or not (i.e., wether the dictionary has been used or not.)
#'            \item \code{level} A numeric value indicating the vertical level of the variable 
#'            (\code{NULL} for 2D variables)
#'      }
#' \item{\code{Data }}{A N-dimensional array. The number of dimensions (N) depends on the 
#' type of request given that dimensions of length one are dropped. Thus, N can take values 
#' from 4 (several members for a rectangular domain with different values for longitude, latitude, 
#' ensemble and time dimensions) to 1 (atomic vector), for single-point and single-member selections, 
#' for which only the time dimension is required. The dimensions are labelled by the \dQuote{dimnames} attribute. 
#' Note that the order of the dimensions is not fixed.}
#' \item{\code{xyCoords }}{A list with \code{x} and \code{y} components, as required by many standard 
#' mapping functions in R (see \code{\link[grDevices]{xy.coords}}. In addition, the attribute \code{projection} 
#' provides geo-referencing information as stored in the original dataset.}
#' \item{\code{Dates }}{A list with two \code{\link[base]{POSIXct}} time elements of the same length as the
#'  \sQuote{time} dimension in \code{Data}, defining the time boundaries of the time axis coordinates in the
#'   interval \emph{[start, end)}, or if the loaded field is static, a character string indicating it.}
#'   
