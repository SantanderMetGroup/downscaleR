#' @importFrom utils packageDescription
#' @importFrom RCurl getURL

.onAttach <- function(...) {
      pkgname <- "downscaleR"
      lib <- system.file(package = pkgname)
      ver <- packageDescription(pkgname)$Version
      builddate <- packageDescription(pkgname)$Date
      mess <- paste(pkgname, " version ", ver, " (", builddate,") is loaded", sep = "")
      packageStartupMessage(mess)
      url <- "https://raw.githubusercontent.com/SantanderMetGroup/downscaleR/stable/DESCRIPTION"
      a <- getURL(url, ssl.verifypeer = FALSE)
      b <- readLines(textConnection(a))
      latest.ver <- package_version(gsub("Version: ", "", b[grep("Version", b)]))
      if (ver < latest.ver) {
            ver.mess1 <- paste("WARNING: Your current version of ", pkgname, " (v", ver, ") is not up-to-date", sep = "")
            ver.mess <- paste("Get the latest stable version (", latest.ver, ") using <devtools::install_github('SantanderMetGroup/downscaleR@stable'>)", sep = "")
            ver.mess
            packageStartupMessage(ver.mess1)
            packageStartupMessage(ver.mess)
      }
}
# End

