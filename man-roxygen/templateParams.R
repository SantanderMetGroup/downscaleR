#' @param dataset A character string indicating the database to be accessed. This is usually a path to a local file or a URL 
#' pointing to a netCDF or NcML file in the case of netCDF and/or gridded datasets. For station data in standard ASCII format,
#' this is the path to the directory the dataset lives in.
#' @param var Variable code (character string). This is the name of the variable according to the R standard naming
#'  (see the next argument). For variables with vertical levels, the vertical level is specified next to the variable name followed
#'   by the \dQuote{@@} symbol (e.g. \code{var = "z@@700"} for geopotential heigth at 700 mb isobaric surface pressure level).
#'   It is also possible to enter the variable name as originally coded in the dataset to skip data homogenization.
#' @param lonLim Vector of length = 2, with minimum and maximum longitude coordinates, in decimal degrees, of the bounding box selected.
#'  For single-point queries, a numeric value with the longitude coordinate. If \code{NULL} (default), the whole longitudinal range
#'   is selected (Note that this may lead to a large output object size).
#' @param latLim Same as \code{lonLim}, but for the selection of the latitudinal range.
#' @param season An integer vector specifying the desired season (in months, January = 1 ..., December = 12).
#'  Options include one to several (contiguous) months. Default to \code{NULL}, indicating a full year selection (same as \code{season = 1:12}).
#' @param years Optional vector of years to select. Default (\code{NULL}) to all available years. If the requested variable is static (e.g. orography)
#'  it will be ignored.  
