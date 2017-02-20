#' @importFrom utils packageDescription
#' @importFrom RCurl getURL

.onAttach <- function(...) {
      pkgname <- "downscaleR"
      lib <- system.file(package = pkgname)
      ver <- packageDescription(pkgname)$Version
      builddate <- packageDescription(pkgname)$Date
      mess <- paste(pkgname, " version ", ver, " (", builddate,") is loaded", sep = "")
      packageStartupMessage(mess)
      url <- paste0("https://raw.githubusercontent.com/SantanderMetGroup/", pkgname, "/master/DESCRIPTION")
      con <- tryCatch(getURL(url, ssl.verifypeer = FALSE), error = function(er) {
            er <- NULL
            return(er)
      })
      # packageStartupMessage("###################################################################\nIMPORTANT NOTE: Since downscaleR v2.0-0, the new dependency package\ntransformeR is introduced. Data manipulation functions are now part\nof transformeR: <https://github.com/SantanderMetGroup/transformeR>\n###################################################################")
      if (!is.null(con)) {
            b <- readLines(textConnection(con))
            latest.ver <- package_version(gsub("Version: ", "", b[grep("Version", b)]))
            if (ver < latest.ver) {
                  ver.mess1 <- paste0("WARNING: Your current version of ", pkgname, " (v", ver, ") is not up-to-date")
                  ver.mess <- paste0("Get the latest stable version (", latest.ver,
                                     ") using <devtools::install_github(c('SantanderMetGroup/transformeR','SantanderMetGroup/", pkgname, "'))>")
                  packageStartupMessage(ver.mess1)
                  packageStartupMessage(ver.mess)
            }
      }
}
# End

