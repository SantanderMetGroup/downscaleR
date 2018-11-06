#' @title Plots daily/annual series and the annual correlation map of different grid objects
#' @description Plots daily/annual series and the annual correlation map of different grid objects
#' @param obs grid of observations. 
#' @param raw grid of model data.
#' @param downscaled Optional. grid of the downscaling output.  
#' @param location Coordinates of a location in the geographic domain of 'obs'. 
#' If NULL it is randomly chosen from coordinates in 'obs'.
#' @param type Character value, either \code{"daily"} or \code{"interannual"}, indicating is the assessment is to
#' be performed on a daily or interannual basis.
#' @param members An integer vector indicating \strong{the position} of the members to be subset. 
#' If NULL all members are considered when type = "interannual" and one member is randomly chosen when type = "daily".
#' @param na.tolerance proportion of NAs in a grid cell (location) that are allowed to calculate correlation. 
#' @param ylim 'ylim' argument passed to the time series plot.
#' @param main 'main' argument passed to the plot.
#' @family visualization
#' @importFrom gridExtra grid.arrange
#' @return Two diagnostic plots with observed, simulated and (possibly) downscaled time series, and a QQ-plot by percentlies.
#' @author M. Iturbide 
#' @export
#' @examples {
#' data("VALUE_Iberia_pr")
#' data("NCEP_Iberia_pr")
#' y <- VALUE_Iberia_pr
#' x <- NCEP_Iberia_pr
#' x$Data <- x$Data*86400
#' quickDiagnostics(obs = y, raw = x, location = c(-2, 43), type = "daily")
#' #quickDiagnostics(obs = y, raw = x, location = c(-2, 43), type = "interannual")
#' eqm1win <- biasCorrection(y = y, x = x,
#'                           method = "eqm",
#'                           precipitation = TRUE,
#'                           extrapolation = "none")
#' quickDiagnostics(y, x, eqm1win)
#' quickDiagnostics(y, x, eqm1win, location = c(-2, 43))
#' #quickDiagnostics(obs = y, raw = x, downscaled = eqm1win,
#' # location = c(-2, 43), type = "interannual")
#' }


quickDiagnostics <- function(obs, 
                             raw, 
                             downscaled = NULL, 
                             location = NULL, 
                             type = c("daily", "interannual"), 
                             members = NULL, 
                             na.tolerance = .3, 
                             ylim = NULL, 
                             main = NULL){
      
      if (is.null(location)) {
            coordref <- getCoordinates(obs)
            if (class(coordref) == "list") coordref <- expand.grid(coordref$x, coordref$y)
            location <- as.numeric(coordref[sample(1:nrow(coordref), 1),])
      }
      type <- match.arg(type, choices = c("daily", "interannual"))
      message("[", Sys.time(), "] Performing quick diagnostics of ", type, " data...")
      if (is.null(members) & "member" %in% names(getShape(raw))) {
            if (type == "daily") {
                  members <- sample(1:getShape(raw)["member"],1)
                  message("[", Sys.time(), "] Member ", members, " is randomly chosen.")
            } else {
                  members <- 1:getShape(raw)["member"]
            }
      }
      suppressWarnings(raw <- subsetGrid(raw, members = members))
      if (!is.null(downscaled)) suppressWarnings(downscaled <- subsetGrid(downscaled, members = members))
      if (type == "daily") {
            if (!is.null(downscaled)) {
                  if (difftime(downscaled$Dates$start[2], downscaled$Dates$start[1], units = "weeks") > 1) {
                        stop("downscaled data is not daily, try with type = 'interannual'")
                  }
            }
            dailyOutlook(obs, raw, downscaled, location, ylim)
      } else if (type == "interannual") {
            interannualOutlook(obs, raw, downscaled, location, na.tolerance = na.tolerance, ylim, main)
      }
      message("[", Sys.time(), "] Done.")
}
#end


#' @importFrom stats cor
#' @importFrom graphics par plot axis text polygon lines legend points
#' @importFrom stats sd cor
#' @importFrom grDevices rgb
#' @keywords internal
#' @importFrom transformeR subsetGrid getCoordinates getYearsAsINDEX draw.world.lines aggregateGrid climatology plotClimatology interpGrid getGrid
#' @importFrom sp SpatialPoints

interannualOutlook <- function(obs, raw, downscaled = NULL, location = c(-42.5, -3), na.tolerance = .3, ylim = NULL, main = NULL){
      period <- unique(getYearsAsINDEX(obs))
      if (!identical(period, unique(getYearsAsINDEX(raw)))) stop("obs and raw have different temporal (year) domain")
      suppressMessages(obsy <- aggregateGrid(obs, aggr.y = list(FUN = "mean", na.rm = T)))
      suppressMessages(rawy <- aggregateGrid(raw, aggr.y = list(FUN = "mean", na.rm = T)))
      if (!is.null(downscaled)) {
            if (!identical(period, unique(getYearsAsINDEX(downscaled)))) stop("obs and downscaled have different temporal (year) domain")
            suppressMessages(downy <- aggregateGrid(downscaled, aggr.y = list(FUN = mean, na.rm = T)))
      }
      if ("member" %in%  getDim(raw)) {
         nmem <- getShape(raw)["member"]
         suppressMessages(rawsd <- aggregateGrid(rawy, aggr.mem = list(FUN = "sd", na.rm = T)))
         suppressMessages(rawy <- aggregateGrid(rawy, aggr.mem = list(FUN = "mean", na.rm = T)))
         if (!is.null(downscaled)){
               suppressMessages(downsd <- aggregateGrid(downy, aggr.mem = list(FUN = "sd", na.rm = T)))
               suppressMessages(downy <- aggregateGrid(downy, aggr.mem = list(FUN = "mean", na.rm = T)))
         } 
      }
      if ("loc" %in% getDim(obs)) {
            a <- abs(getCoordinates(obs)$y - location[2]) + abs(getCoordinates(obs)$x - location[1])
            indstation <- which(a == min(a))
            x <- obsy$Data[,indstation]
            location[1] <- obs$xyCoords[indstation,1]
            location[2] <- obs$xyCoords[indstation,2]
            warning("Coordinates of the nearest station to the specified location are considered: 
                    x = ", location[1], ", y = ", location[2])
            if (!is.null(downscaled)) {
                  w <- downy$Data[,indstation]
                  if ("member" %in%  getDim(downscaled)) wsd <- downsd$Data[,indstation]
            } 
      } else {
            x <- subsetGrid(obsy, lonLim = location[1], latLim = location[2])$Data
            if (!is.null(downscaled)) {
                  w <- subsetGrid(downy, lonLim = location[1], latLim = location[2])$Data
                  if ("member" %in%  getDim(downscaled)) wsd <- subsetGrid(downsd, lonLim = location[1], latLim = location[2])$Data
            }
      }
      y <- subsetGrid(rawy, lonLim = location[1], latLim = location[2])$Data
      if ("member" %in%  getDim(raw)) ysd <- subsetGrid(rawsd, lonLim = location[1], latLim = location[2])$Data
      #statistics
      rmse.direct <- sqrt(mean((x - y)^2, na.rm = T))
      bias.direct <-  sum(y - x, na.rm = T)/sum(x,  na.rm = T)
      rho.direct <- cor(x = x, y = y, method = "spearman")
      if (!is.null(downscaled)){
            rmse.down <- sqrt(mean((x - w)^2, na.rm = T))
            bias.down <-  sum(w - x, na.rm = T)/sum(x, na.rm = T)
            rho.down <- cor(x = x, y = w, method = "spearman")
      }
      if (is.null(ylim)) {
            mi <- floor(min(c(x, y))) - 1
            ma <-  floor(max(c(x, y)))
            ylim <- c(mi, ma + (ma - mi))
      }
      plot(1:length(period), x, xlim = c(0,length(period)), ylim = ylim, xlab="", xaxt = "n", 
           ylab = "Annual/seasonal mean value", cex = .6, col = NULL, main = main)
      tck <- axis(1, at = 1:length(period), labels=FALSE)
      text(tck,  par("usr")[3] - 2, xpd = TRUE, labels = period, 
           srt = 90, cex =.6)
      #plot the sd (shadows)
      lines(1:length(period), x, col = "black")
      lines(1:length(period), y, col = "red")
      if (!is.null(downscaled)) lines(1:length(period), w, col = "blue")
      if ("member" %in%  getDim(raw)) {
            polygon(x = c(1:length(period), length(period):1), y = c(ysd + y,rev(y - ysd)), col = rgb(1,0,0,0.2), border = NA)
            if (!is.null(downscaled)) {
                  polygon(x = c(1:length(period), length(period):1), y =c (wsd+w,rev(w-wsd)), col = rgb(0,0,1,0.2), border = NA)
            }
      }
      #legend
      if (!is.null(downscaled)) {
            legend(0, ylim[2] , legend = c("obs",  
                                           paste("direct: rho=", 
                                                 round(rho.direct, digits = 2), ", bias=", 
                                                 round(bias.direct, digits = 2), sep = ""), 
                                           paste("downscaled: rho=", 
                                                 round(rho.down, digits = 2), ", bias=", 
                                                 round(bias.down, digits = 2), sep = "")), 
                         fill = c("black", "red", "blue"), box.lwd = 0, cex = .8)
       } else {
            legend(0, ylim[2] , legend = c("obs",  
                                           paste("direct: rho=", 
                                                 round(rho.direct, digits = 2), ", bias=", 
                                                 round(bias.direct, digits = 2), sep = "")), 
                   fill = c("black", "red"), box.lwd = 0, cex = .8)
       }
       pl <- grabGrob()
       ## correlation map
       suppressMessages(corgrid <- climatology(obs))
       if (is.null(downscaled)) {
             downy <- suppressMessages(interpGrid(rawy, getGrid(obsy)))
       }
       if ("loc" %in% getDim(obs)) {
             corgrid$Data[1,] <- vapply(1:nrow(getCoordinates(corgrid)), FUN.VALUE = numeric(length = 1), FUN = function(m) {
                   cor(obsy$Data[,m], downy$Data[,m], method = "spearman") 
             })
       } else {
             for (i in 1:length(getCoordinates(corgrid)$x)) {
                   corgrid$Data[1,,i] <- vapply(1:length(getCoordinates(corgrid)$y), FUN.VALUE = numeric(length = 1), FUN = function(m) {
                         cor(obsy$Data[,m,i], downy$Data[,m,i], method = "spearman")  
                  })
             }
       }
       plcor <- 
             plotClimatology(corgrid, backdrop.theme = "countries", cuts = seq(-1,1,0.25),
                                key.space = "bottom", scales=list(draw = TRUE),
                             auto.key = T,par.settings=list(fontsize=list(text= 8)),
                       sp.layout = list(list(SpatialPoints(matrix(location,ncol = 2)), 
                                             first = FALSE, pch = 2, cex = 1.8, col = "black")))
       grid.arrange(pl, plcor, nrow = 1, ncol = 2, widths = c(1.4,1))
}
#end


#' @importFrom stats cor quantile
#' @importFrom graphics par plot lines legend points
#' @importFrom stats sd cor
#' @importFrom grDevices rgb
#' @importFrom transformeR subsetGrid
#' @keywords internal

dailyOutlook <- function(obs, raw, downscaled = NULL, location = c(-42.5, -3), ylim = NULL){
      if ("loc" %in% getDim(obs)) {
            a <- abs(getCoordinates(obs)$y - location[2]) + abs(getCoordinates(obs)$x - location[1])            
            indstation <- which(a == min(a))
            x <- obs$Data[,indstation]
            location[1] <- obs$xyCoords[indstation,1]
            location[2] <- obs$xyCoords[indstation,2]
            warning("Coordinates of the nearest station to the specified location are considered: 
                    x = ", location[1], ", y = ", location[2])
      } else {
            x <- subsetGrid(obs, lonLim = location[1], latLim = location[2], outside = T)$Data
      }
      y <- subsetGrid(raw, lonLim = location[1], latLim = location[2], outside = T)$Data
      if (!is.null(downscaled)) {
            if ("loc" %in% getDim(downscaled)) {
                  w <- downscaled$Data[,indstation]      
            } else {
                  w <- subsetGrid(downscaled, lonLim = location[1], latLim = location[2], outside = T)$Data
            }
      }
      yran <- if (is.null(ylim)) {
            mi <- 0
            ma <-  floor(max(c(x,y), na.rm = T))
            c(mi, (ma + ma/4))
      } else {
            ylim
      }
      # Daily time series plot
      par(mfrow=c(1,2))
      plot(1:length(x), 
           x, 
           lwd = 2, ty = "l", lty = 3,  
           ylim = yran, xlim = c(1,dim(obs$Data)[1]),
           xlab = "Days", ylab = obs$Variable$varName, main = "Daily series")
      rmse.direct <- sqrt(mean((x - y)^2, na.rm = T))
      bias.direct <-  sum(y - x, na.rm = T) / sum(x, na.rm = T)
      naind <- intersect(which(!is.na(x)), which(!is.na(y)))
      rho.direct <- cor(x = x[naind], y = y[naind], method = "spearman")
      lines(1:length(y), y, col = "red", lwd = 1)
      if (!is.null(downscaled)) {
            xt <- x[(length(x) - (length(w) - 1)):length(x)]
            yt <- y[(length(x) - (length(w) - 1)):length(x)]
            rmse.down <- sqrt(mean((xt - w)^2, na.rm = T))
            bias.down <-  sum(w - xt, na.rm = T)/sum(x, na.rm = T)
            naind2 <- intersect(which(!is.na(xt)), which(!is.na(yt)))
            rho.down <- cor(x = xt[naind2], y = w[naind2], method = "spearman")
            rmse.direct <- sqrt(mean((xt - yt)^2, na.rm = T))
            bias.direct <-  sum(yt - xt, na.rm = T)/sum(x, na.rm = T)
            rho.direct <- cor(x = xt[naind2], y = yt[naind2], method = "spearman")
            lines((length(x) - (length(w) - 1)):length(x), w, col = "blue", lwd = 1)
                  legend(0, yran[2] , legend = c("obs",  
                                                 paste("raw: rho=", 
                                                       round(rho.direct, digits = 2), ", bias=", 
                                                       round(bias.direct, digits = 2), sep = ""), 
                                                 paste("downscaled: rho=", 
                                                       round(rho.down, digits = 2), ", bias=", 
                                                       round(bias.down, digits = 2), sep = "")), 
                         fill = c("black", "red", "blue"), box.lwd = 0, cex = .6)
       } else {
            legend(0, yran[2] , legend = c("obs", paste("raw: rho=", 
                                            round(rho.direct, digits = 2), ", bias=", 
                                            round(bias.direct, digits = 2), sep = "")), 
                         fill = c("black", "red"), box.lwd = 0, cex = .6)
      }
      # qq-plot
      qy <- quantile(y, probs = seq(0.01, .99, 0.01), na.rm = T, , type =4)
      if (!is.null(downscaled)){
            qw <-  quantile(w, probs = seq(0.01, .99, 0.01), na.rm = T, type =4)
      }else{
            qw <- NA
      }
      q1 <- quantile(x, probs = seq(0.01, .99, 0.01), na.rm = T, , type =4)
      yran <- c(0, max(c(q1, qy, qw), na.rm = T))
      plot(q1, qy, col="red", main = "qq-plot", xlab = "obs", ylab = "predicted", ylim = yran)
      lines(0:max(c(q1, qy, qw), na.rm = T), 0:max(c(q1, qy, qw), na.rm = T))
      if (!is.null(downscaled)){
            points(q1,qw, col="blue")
      }
      par(mfrow = c(1,1)) 
}

#end    


grabGrob <- function(){
      grid.echo()
      grid.grab()
}