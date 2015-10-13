#' @title Login to the Santander Met Group Thredds Data Server
#' @description Provides HTTP authentication for accessing the Santander Met Group's THREDDS Data Server
#' @param username A character string with a valid user ID. See details.
#' @param password Character string. Authorized password. See details.
#' @param proxy.host In case of a proxy connection, the name of the proxy host, as a character string
#' @param proxy.port In case of a proxy connection, an integer indicating the proxy port
#' @details The Santander Met Group has deployed a THREDDS Data Server to enable access to different climate databases/data collections,
#' implementing a fine-grained user authorization and using remote data access protocols with data subsetting capabilities.
#' Prior to data access, users must log in. Registration can be obtained via the THREDDS Administration Panel (\href{http://www.meteo.unican.es/tap}{TAP}),
#'  indicating the group (Project) you belong to (e.g. CORDEX, SPECS ...), which will grant access to certain databases.
#'  Further details on registration for data access can be viewed in this \href{http://meteo.unican.es/ecoms-udg/DataServer/Registration}{example link}.
#' @author J Bedia \email{joaquin.bedia@@gmail.com}, M. Vega.
#' @export


loginUDG <- function(username, password, proxy.host = NULL, proxy.port = NULL) {
      proxy.port <- as.integer(proxy.port)
      if (!is.character(username) | !is.character(password)) {
            stop("\'username\' and \'password\' must be character strings")
      }
      if (!is.null(proxy.host)) {
            J("ucar.nc2.util.net.HTTPSession")$setGlobalProxy(proxy.host, proxy.port)
      }
      J("ucar.httpservices.MyHTTPFactory")$setCredentials(username, password)
}
# End

