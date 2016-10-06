
##  It produces a  However, the original ecomsUDG data is saved on a 6h time step which is aggregated into daily before creating the field object. This is mainly due to the fact that the class of the time index is usually Date. So it does not deal with subdaily data. One way to overcome this limitation is to work with POSIXct and POSIXlt classes to represent calendar dates and times. (see also ? strptime) 
#' @title Conversion to esd (package) object
#' @description The udg2esd function is a conversion tool between the ecomsUDG server and the ESD package.
#' @param x grid  
#' @return grid object that can be further analysed using the {esd} functionalities.
#' @author M. Iturbide 
#' @export


udg2esd <- function(x,verbose=FALSE) {
      ## load ecomUDG.Raccess library
      ##
      udg3d2field <- function(x) {
            ## Permute dimensions to lon,lat,time
            xp <- aperm(x,c(3,2,1))
            ##xip <- xi
            d <- dim(xp)
            dim(xp) <- c(d[1]*d[2],d[3])
            z <- zoo(t(xp),order.by=index)
            y <- as.field(x=z,unit=unit,lon=lon,lat=lat,param=param,
                          longname=longname,src =src,info=info)
            invisible(y)
      }
      param <- x$Variable$varName
      unit <- attr(x$Variable, "units")
      longname <- attr(x$Variable, "longname")
      src <- attributes(x)$dataset
      info <- attr(x,'source')
      ## browser()
      ## convert sp coordinates to lon and lat
      lonlat <- as.data.frame(x$xyCoords)
      lon <- lonlat$x
      lat <- lonlat$y
      ## Get time index
      index <- as.Date(x$Dates$start)
      ## Get dimensions
      d <- dim(x$Data)
      ##
      if (!is.null(dim(x))) n <- dim(x)[2] else n <- 1
      ## browser()
      ##for (i in 1 : n) { ## change this value later AM 05-05-2014
      ## convert to zoo object
      xi <- x$Data
      n <- length(d)
      if (verbose)
            print(paste('Dimensions:',paste(dim(xi),collapse='/')),sep='')
      ## aggregate into daily values
      ## z <- aggregate(z,FUN=mean,by=as.Date(index))
      ## print("Aggregate to daily values")
      if (n==0) {## convert as station object if single grid value
            if ((length(lon)==1) & (length(lat)==1)) {
                  z <- zoo(xi,order.by=index)
                  y <- as.station(x=z,loc= names(x$MemberData)[i],
                                  unit=unit,lon=lon,lat=lat,param=param,
                                  longname=longname,
                                  src =attributes(x)$dataset,
                                  info=attr(x,'source'))
                  if (i == 1) Y <- y else {browser() ; Y <- combine(y,Y)}
            }
      }
      else if (n==3) { ## Convert as field object
            y <- udg3d2field(x$Data)
      } else if (n==4) {
            mem <- grep('memb',attr(x$Data,"dimensions"))
            n.len <- d[mem]
            y <- lapply(1:n.len, function(i) {
                  mem <- x$Data[i,,,]
                  udg3d2field(mem)
            })##y <- udg4d2field(x)
            names(y) <- x$Members
      } else stop('x has dimension longer than 4, do not know how to handle')
      invisible(y)
}