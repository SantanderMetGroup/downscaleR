##########################3
#' Java string dates to R character
#' 
#' Converts objects of the Java class \sQuote{java/lang/String} to character in R
#' 
#' @param javaString A \sQuote{java/lang/String} object or array 
#' @return A character vector in R
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}

javaString2rChar <- function(javaString) {
      r.string <- unlist(strsplit(gsub("\\[|]|\\s", "", javaString), ","))
      return(r.string)
}
# End