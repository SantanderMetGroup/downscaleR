#' @section Geolocation parameters: Regarding the selection of the spatial domain,
#'  it is possible to select the whole spatial domain of the dataset by defining the arguments \code{lonLim=NULL}
#'  and \code{latLim=NULL}. More often, rectangular domains are defined by the minimum and maximum coordinates
#'  in longitude and latitude (for instance \code{lonLim=c(-10,10)} and \code{latLim=c(35,45)} indicates a
#'  rectangular window centered in the Iberian Peninsula), and single grid-cell values
#'  (for instance \code{lonLim=-3.21} and \code{latLim=41.087} for retrieving the data in the closest grid
#'  point to the point coordinate -3.21E, 41.087N. In the last two cases, the function
#'  operates by finding the nearest (euclidean distance) grid-points to the coordinates introduced.
#'  
#'  In the case of station data (\code{\link{loadStationData}}), the logic is the same, taking into account that in the case
#'  of rectangular domains, all stations falling inside that window will be loaded. For single-point selections,
#'  the closest station will be chosen, and a note on-screen will inform about the distance from the selected point
#'  to the chosen station.
#'  
