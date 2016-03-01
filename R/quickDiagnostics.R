#' @title Plots daily/annual series and the annual correlation map of different grid objects
#' @description Plots daily/annual series and the annual correlation map of different grid objects
#' @param obs grid of observations. 
#' @param sim grid of model data.
#' @param downscaled Optional. grid of the downscaling output.  
#' @param location Coordinates of a location in the geographic domain of the grid.
#' @param na.tolerance proportion of NAs in a grid cell (location) that are allowed to calculate correlation. 
#' @param ylim 'ylim' argument passed to the time series plot.
#' @param main 'main' argument passed to the plot.
#' @param type Character value, either \code{"daily"} or \code{"interannual"}, indicating is the assessment is to
#' be performed on a daily or interannual basis.
#' @family visualization
#' @importFrom fields image.plot
#' @return Two diagnostic plots with observed, simulated and (possibly) downscaled time series, and a QQ-plot by percentlies.
#' @author M. Iturbide 
#' @export


quickDiagnostics <- function(obs, sim, downscaled = NULL, location = c(-42.5, -3), type = c("daily", "interannual"), na.tolerance = .3, ylim = NULL, main = NULL){
      
      if (type == "daily") {
            if (!is.null(downscaled)) {
                  if(difftime(downscaled$Dates$start[2], downscaled$Dates$start[1], units = "weeks") > 1){
                        stop("downscaled data is not daily, try with type = 'interannual'")
                  }
            }
            dailyOutlook (obs, sim, downscaled, location, ylim)
      } else if (type == "interannual") {
            interannualOutlook (obs, sim, downscaled, location, na.tolerance = na.tolerance, ylim, main)
      }
}
#end

interannualOutlook <- function(obs, sim, downscaled = NULL, location = c(-42.5, -3), na.tolerance = .3, ylim = NULL, main = NULL){
      par(mfrow = c(1,2))
      period.id <- (getYearsAsINDEX(sim))
      period <- unique(period.id)
      if (!is.null(downscaled)){
            test.id2 <- getYearsAsINDEX(downscaled)
            train.id <- period.id[which(is.na(match(period.id, test.id2)))]
            if(length(train.id)==0){train.id <-period.id}
            test.id <- period.id[which(is.na(match(period.id, train.id)))]
            if(length(test.id)==0){test.id <- period.id}
            test.period <- unique(test.id)
            #             train.period <- period[which(is.na(match(period,test.period)))]
            train.period <-unique(train.id)
            obs.test <- subsetGrid(obs, years = test.period) 
            obs.train <- subsetGrid(obs, years = train.period)
            sim.test <- subsetGrid(sim, years = test.period)
            sim.train <- subsetGrid(sim, years = train.period)
            if(length(which(is.na(match(train.id, test.id2))))==0){
                  comper <- TRUE
            }else{comper <- FALSE}
      }
      x.coord <- getCoordinates(obs)$x
      y.coord <- getCoordinates(obs)$y
      coords <- expand.grid(1:length(x.coord), 1:length(y.coord))
      xo <- findInterval(location[1], x.coord)
      yo <- findInterval(location[2], y.coord)
      xi <- findInterval(x.coord[xo], getCoordinates(sim)$x)
      yi <- findInterval(y.coord[yo], getCoordinates(sim)$y)
      if (any(attr(sim$Data, "dimensions")=="member")){
            nmem <- dim(sim$Data)[which(attr(sim$Data, "dimensions")=="member")]
            #Daily time series plot
            ##Compute statistic
            x <- tapply(obs$Data[, yo,xo], INDEX = period.id, FUN = mean)
            yl <- lapply(1:nmem, function(m) {
                  tapply(sim$Data[m,,yi,xi], INDEX = period.id, FUN = mean)
            })
            y <- Reduce("+", yl)/length(yl)
            ys <- apply(simplify2array(yl), 1, sd)
            if(!is.null(downscaled)){
                  xt <- tapply(obs.test$Data[, yo,xo], INDEX = test.id, FUN = mean)
                  ytl <- lapply(1:nmem, function(m){
                        tapply(sim.test$Data[m,,yi,xi], INDEX = test.id, FUN = mean)})
                  yt <- Reduce("+", ytl)/length(ytl)
                  wl <- lapply(1:nmem, function(m){
                        if(difftime(downscaled$Dates$start[2], downscaled$Dates$start[1], units = "weeks") > 1){
                              downscaled$Data[m,,yo,xo]   
                        }else{
                              tapply(downscaled$Data[m,,yo,xo], INDEX = test.id2, FUN = mean)}
                  })
                  
                  w <- Reduce("+", wl)/length(wl)
                  ws <- apply(simplify2array(wl), 1, sd)
                  ## ACCURACY
                  ### Spearman correlation rho and the Root Mean Square Error (RMSE) as accuracy measures 
                  ### for the direct and calibrated simulation in the TEST PERIOD
                  rmse.down <- sqrt(mean((xt - w)^2))
                  bias.down <-  sum(w - xt)/sum(x)
                  rho.down <- cor(x = xt, y = w, method = "spearman")
                  rmse.direct <- sqrt(mean((xt - yt)^2))
                  bias.direct <-  sum(yt - xt)/sum(x)
                  rho.direct <- cor(x = xt, y = yt, method = "spearman")
            } else {
                  rmse.direct <- sqrt(mean((x - y)^2))
                  bias.direct <-  sum(y - x)/sum(x)
                  rho.direct <- cor(x = x, y = y, method = "spearman")
            }
            if (is.null(ylim)) {
                  mi <- floor(min(c(x,y)))-1
                  ma <-  floor(max(c(x,y)))
                  ylim <- c(mi, ma + (ma-mi))
            }
            plot(1:length(period), x, xlim = c(0,length(period)), ylim = ylim, xlab="", xaxt = "n", 
                 ylab = "Annual/seasonal mean value", cex = .6, col = NULL, main = main)
            tck <- axis(1, at = 1:length(period), labels=FALSE)
            text(tck,  par("usr")[3] - 2, xpd = TRUE, labels = (1981:2010), 
                 srt = 90, cex =.6)
            #plot the sd (shadows)
            polygon(x = c(1:length(period), length(period):1), y =c (ys+y,rev(y-ys)), col = rgb(1,0,0,0.2), border = NA)
            if (!is.null(downscaled)) {
                  if(comper == TRUE){
                        polygon(x = c(1:length(period), length(period):1), 
                                y =c (ws+w,rev(w-ws)), col = rgb(0,0,1,0.2), border = NA)
                        #plot the mean (lines)
                        lines(1:length(period), x[1:length(period)], lwd = 2, xlim = c(0,length(period)))
                        lines(1:length(period), y[1:length(period)], lwd = 2, xlim = c(0,length(period)),
                              col="red")
                        lines(1:length(period), w, col = "blue", lwd = 2)     
                  }else{
                        polygon(x = c((length(train.period)+1):length(period), length(period):(length(train.period)+1)), 
                                y =c (ws+w,rev(w-ws)), col = rgb(0,0,1,0.2), border = NA, main = main)
                        #plot the mean (lines)
                        lines(1:(length(train.period)+1), x[1:(length(train.period)+1)], lwd = 2, lty = 4, xlim = c(0,length(period)))
                        lines((1+length(train.period)):length(period) , x[(1+length(train.period)):length(period)], 
                              lwd = 2)
                        lines(1:(length(train.period)+1), y[1:(length(train.period)+1)], lwd = 2, lty = 4, xlim = c(0,length(period)),
                              col="red")
                        lines((1+length(train.period)):length(period) , y[(1+length(train.period)):length(period)], 
                              lwd = 2, col ="red")
                        lines((length(train.period)+1):length(period), w, col = "blue", lwd = 2)
                  }
                  legend(0, ylim[2] , legend = c("obs",  
                                                 paste("direct: rho=", 
                                                       round(rho.direct, digits = 2), ", bias=", 
                                                       round(bias.direct, digits = 2), sep = ""), 
                                                 paste("downscaled: rho=", 
                                                       round(rho.down, digits = 2), ", bias=", 
                                                       round(bias.down, digits = 2), sep = "")), 
                         fill = c("black", "red", "blue"), box.lwd = 0, cex = .8)
            } else {
                  # plot the mean (lines)
                  lines(1:length(period), x, lwd = 2)
                  lines(1:length(period), y, lwd = 2, col ="red")
                  legend(0, ylim[2] , legend = c("obs",  
                                                 paste("direct: rho=", 
                                                       round(rho.direct, digits = 2), ", bias=", 
                                                       round(bias.direct, digits = 2), sep = "")), 
                         fill = c("black", "red"), box.lwd = 0, cex = .8)
            }
            ## correlation map
            r <- matrix(nrow = length(x.coord), ncol = length(y.coord))
            for (i in 1:nrow(coords)) {
                  if(!is.null(downscaled)){
                        x <- tapply(obs.test$Data[, coords[i,2], coords[i,1]],
                                    INDEX = test.id, FUN = mean)
                        
                        wl <- lapply(1:nmem, function(m) {
                              if(difftime(downscaled$Dates$start[2], downscaled$Dates$start[1], units = "weeks") > 1){
                                    downscaled$Data[m,,coords[i,2], coords[i,1]]   
                              }else{
                                    tapply(downscaled$Data[m,, coords[i,2], coords[i,1]],
                                           INDEX = test.id, FUN = mean)}
                        })
                  } else {
                        x <- tapply(obs$Data[, coords[i,2], coords[i,1]],
                                    INDEX = period.id, FUN = mean)
                        wl <- lapply(1:nmem, function(m) {
                              
                              tapply(sim$Data[m,, coords[i,2], coords[i,1]],
                                     INDEX = period.id, FUN = mean)
                        }) 
                  }
                  w <- Reduce("+", wl)/length(wl)
                  r[coords[i,1],coords[i,2]] <- cor(x,w,method = "spearman")
            }
            # Ignore negative values
            r[which(r < 0)] <- 0
            image.plot(x.coord, y.coord, r, asp = 1, breaks = seq(0,1,0.1),
                       nlevel = 10, lab.breaks = c("=<0",as.character(seq(0.1,1,0.1))),
                       xlab = "longitude", ylab = "latitude")
            
            
            
            points(location[1], location[2], cex = 2, pch = 17)
      } else {
            ##Compute statistic
            x <- tapply(obs$Data[, yo, xo], INDEX = period.id, FUN = mean)
            y <- tapply(sim$Data[, yi, xi], INDEX = period.id, FUN = mean)
            
            if (!is.null(downscaled)) {
                  
                  if(difftime(downscaled$Dates$start[2], downscaled$Dates$start[1], units = "weeks") > 1){
                        w <- downscaled$Data[,yo,xo]   
                  }else{
                        w <- tapply(downscaled$Data[,yo,xo], INDEX = test.id2, FUN = mean)
                        
                  }
                  if(comper == TRUE){
                        xt <- x
                        yt <- y
                  }else{
                        xt <- x[(1+length(train.period)):length(period)]
                        yt <- y[(1+length(train.period)):length(period)]
                  }
                  ## ACCURACY
                  ### Spearman correlation rho and the Root Mean Square Error (RMSE) as accuracy measures 
                  ### for the direct and calibrated simulation in the TEST PERIOD
                  if(length(which(is.na(w))) > 0 & (length(which(is.na(w)))/length(w)) < na.tolerance){  
                        ind <- which(is.na(w))
                        rmse.down <- sqrt(mean((xt[-ind] - w[-ind])^2))
                        bias.down <-  sum(w[-ind] - xt[-ind])/sum(x[-ind])
                        rho.down <- cor(x = xt[-ind], y = w[-ind], method = "spearman")
                        rmse.direct <- sqrt(mean((xt - yt)^2))
                        bias.direct <-  sum(yt - xt)/sum(x)
                        rho.direct <- cor(x = xt, y = yt, method = "spearman")
                  }else if(length(which(is.na(w))) == 0){
                        rmse.down <- sqrt(mean((xt - w)^2))
                        bias.down <-  sum(w - xt)/sum(x)
                        rho.down <- cor(x = xt, y = w, method = "spearman")
                        rmse.direct <- sqrt(mean((xt - yt)^2))
                        bias.direct <-  sum(yt - xt)/sum(x)
                        rho.direct <- cor(x = xt, y = yt, method = "spearman")
                  }else{                  
                        stop("Too many NAs in the selected location. Select a new location")
                  }
                  
            } else {
                  rmse.direct <- sqrt(mean((x - y)^2))
                  bias.direct <-  sum(y - x)/sum(x)
                  rho.direct <- cor(x = x, y = y, method = "spearman")
            }
            if (is.null(ylim)) {
                  mi <- floor(min(c(x,y)))-1
                  ma <-  floor(max(c(x,y)))
                  ylim <- c(mi, ma + (ma-mi))}
            plot(1:length(period), x, xlim = c(0,length(period)), ylim = ylim, xlab="", xaxt = "n", 
                 ylab = "Annual/seasonal mean value", cex = .6, col = NULL, main = main)
            tck <- axis(1, at = 1:length(period), labels=FALSE)
            text(tck,  par("usr")[3] - 2, xpd = TRUE, labels = (1981:2010), 
                 srt = 90, cex =.6)
            if (!is.null(downscaled)) {
                  #plot the mean (lines)
                  if(comper == TRUE){
                        lines(1:length(period), x, lwd = 2, xlim = c(0,length(period)))
                        
                        lines(1:length(period), y, lwd = 2, xlim = c(0,length(period)),
                              col="red")
                        lines(1:length(period), w, col = "blue", lwd = 2)                        
                  }else{
                        lines(1:(length(train.period)+1), x[1:(length(train.period)+1)], lwd = 2, lty = 4, xlim = c(0,length(period)))
                        lines((1+length(train.period)):length(period) , x[(1+length(train.period)):length(period)], 
                              lwd = 2)
                        lines(1:(length(train.period)+1), y[1:(length(train.period)+1)], lwd = 2, lty = 4, xlim = c(0,length(period)),
                              col="red")
                        lines((1+length(train.period)):length(period) , y[(1+length(train.period)):length(period)], 
                              lwd = 2, col ="red")
                        lines((length(train.period)+1):length(period), w, col = "blue", lwd = 2)
                  }
                  legend(0, ylim[2] , legend = c("obs",  
                                                 paste("direct: rho=", 
                                                       round(rho.direct, digits = 2), ", bias=", 
                                                       round(bias.direct, digits = 2), sep = ""), 
                                                 paste("downscaled: rho=", 
                                                       round(rho.down, digits = 2), ", bias=", 
                                                       round(bias.down, digits = 2), sep = "")), 
                         fill = c("black", "red", "blue"), box.lwd = 0, cex = .8)
            } else {
                  #plot the mean (lines)
                  lines(1:length(period), x, lwd = 2)
                  lines(1:length(period), y, lwd = 2, col ="red")
                  legend(0, ylim[2] , legend = c("obs",  
                                                 paste("direct: rho=", 
                                                       round(rho.direct, digits = 2), ", bias=", 
                                                       round(bias.direct, digits = 2), sep = "")), 
                         fill = c("black", "red", "blue"), box.lwd = 0, cex = .8)
            }
            #correlation map
            r <- matrix(nrow = length(x.coord), ncol = length(y.coord))
            for (i in 1:nrow(coords)) {
                  if(!is.null(downscaled)) {
                        x <- tapply(obs.test$Data[, coords[i,2], coords[i,1]],
                                    INDEX = test.id, FUN = mean)
                        if(difftime(downscaled$Dates$start[2], downscaled$Dates$start[1], units = "weeks") > 1){
                              w <- downscaled$Data[,coords[i,2],coords[i,1]]   
                        }else{
                              
                              w <- tapply(downscaled$Data[, coords[i,2], coords[i,1]],
                                          INDEX = test.id, FUN = mean)
                        }
                  } else {
                        x <- tapply(obs$Data[, coords[i,2], coords[i,1]],
                                    INDEX = period.id, FUN = mean)
                        w <- tapply(sim$Data[, coords[i,2], coords[i,1]],
                                    INDEX = period.id, FUN = mean)   
                  }
                  if(length(which(is.na(w))) > 0 & (length(which(is.na(w)))/length(w)) < na.tolerance){
                        ind <- which(is.na(w))
                        r[coords[i,1],coords[i,2]] <- cor(x[-ind],w[-ind],method = "spearman")
                        
                  }else{
                        r[coords[i,1],coords[i,2]] <- cor(x,w,method = "spearman")
                        
                  }
                  
            }
            # Ignore negative values
            r[which(r < 0)] <- 0
            image.plot(x.coord, y.coord, r, asp = 1, breaks = seq(0,1,0.1),
                       nlevel = 10, lab.breaks = c("=<0",as.character(seq(0.1,1,0.1))),
                       xlab = "longitude", ylab = "latitude")
            draw.world.lines()
            points(location[1], location[2], cex = 2, pch = 17)
      }
      par(mfrow = c(1,1)) 
}
#end

dailyOutlook <- function(obs, sim, downscaled = NULL, location = c(-42.5, -3), member = NULL, ylim = NULL){
      if(any(attr(obs$Data, "dimensions")=="station")){
            x <- obs$Data
      }else{
            x <- subsetGrid(obs, lonLim = location[1], latLim = location[2])$Data
      }
      y <- subsetGrid(sim, lonLim = location[1], latLim = location[2])$Data
      if(!is.null(downscaled)){
            w <-  subsetGrid(downscaled, lonLim = location[1], latLim = location[2])$Data
      }
      yran <- if (is.null(ylim)) {
            mi <- 0
            ma <-  floor(max(x))
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
           xlab = "Days", ylab = "tp", main = "Daily series")
      
      o <- match(attr(y, "dimensions"),"member")
      mi <- which(!is.na(o))
      
      if(length(mi)>0){
            nmem <- dim(sim$Data)[which(attr(sim$Data, "dimensions")=="member")]
            
            yl <- lapply(1:nmem, function(m){
                  y[m,]
            })
            y <- Reduce("+", yl)/length(yl)
            
            rmse.direct <- sqrt(mean((x - y)^2))
            bias.direct <-  sum(y - x)/sum(x)
            rho.direct <- cor(x = x, y = y, method = "spearman")
            
            lines(1:length(y), 
                  y, 
                  col = "red", lwd = 1)
            if (!is.null(downscaled)){
                  wl <- lapply(1:nmem, function(m){
                        w[m,]
                  })
                  w <- Reduce("+", wl)/length(wl)
                  
                  xt <- x[(length(x) -(length(w)-1)):length(x)]
                  yt <- y[(length(x) -(length(w)-1)):length(x)]
                  rmse.down <- sqrt(mean((xt - w)^2))
                  bias.down <-  sum(w - xt)/sum(x)
                  rho.down <- cor(x = xt, y = w, method = "spearman")
                  rmse.direct <- sqrt(mean((xt - yt)^2))
                  bias.direct <-  sum(yt - xt)/sum(x)
                  rho.direct <- cor(x = xt, y = yt, method = "spearman")
                  lines((length(x) -(length(w)-1)):length(x), 
                        w, col = "blue", lwd = 1)
                  legend(0, yran[2] , legend = c("obs",  
                                                 paste("sim: rho=", 
                                                       round(rho.direct, digits = 2), ", bias=", 
                                                       round(bias.direct, digits = 2), sep = ""), 
                                                 paste("downscaled: rho=", 
                                                       round(rho.down, digits = 2), ", bias=", 
                                                       round(bias.down, digits = 2), sep = "")), 
                         fill = c("black", "red", "blue"), box.lwd = 0, cex = .8)
            } else {
                  legend(0, yran[2] , legend = c("obs",  
                                                 paste("sim: rho=", 
                                                       round(rho.direct, digits = 2), ", bias=", 
                                                       round(bias.direct, digits = 2), sep = "")), 
                         fill = c("black", "red"), box.lwd = 0, cex = .8)
            }
      } else {
            rmse.direct <- sqrt(mean((x - y)^2))
            bias.direct <-  sum(y- x)/sum(x)
            rho.direct <- cor(x = x, y = y, method = "spearman")
            lines(1:length(y), 
                  y, 
                  col = "red", lwd = 1)
            if (!is.null(downscaled)) {
                  xt <- x[(length(x) -(length(w)-1)):length(x)]
                  yt <- y[(length(x) -(length(w)-1)):length(x)]
                  rmse.down <- sqrt(mean((xt - w)^2))
                  bias.down <-  sum(w - xt)/sum(x)
                  rho.down <- cor(x = xt, y = w, method = "spearman")
                  rmse.direct <- sqrt(mean((xt - yt)^2))
                  bias.direct <-  sum(yt- xt)/sum(x)
                  rho.direct <- cor(x = xt, y = yt, method = "spearman")
                  lines((length(x) -(length(w)-1)):length(x), 
                        w, 
                        col = "blue", lwd = 1)
                  legend(0, yran[2] , legend = c("obs",  
                                                 paste("sim: rho=", 
                                                       round(rho.direct, digits = 2), ", bias=", 
                                                       round(bias.direct, digits = 2), sep = ""), 
                                                 paste("downscaled: rho=", 
                                                       round(rho.down, digits = 2), ", bias=", 
                                                       round(bias.down, digits = 2), sep = "")), 
                         fill = c("black", "red", "blue"), box.lwd = 0, cex = .8)
            } else {
                  legend(0, yran[2] , legend = c("obs",  
                                                 paste("sim: rho=", 
                                                       round(rho.direct, digits = 2), ", bias=", 
                                                       round(bias.direct, digits = 2), sep = "")), 
                         fill = c("black", "red"), box.lwd = 0, cex = .8)
            }
      }
      # qq-plot
      q1 <- quantile(x, probs = seq(0.01, .99, 0.01), na.rm = T, , type =4)
      yran <- c(0, max(q1))
      plot(q1, 
           quantile(y, probs = seq(0.01, .99, 0.01), na.rm = T, type =4), 
           col="red", main = "qq-plot", xlab = "obs", ylab = "predicted", ylim = yran)
      lines(0:max(q1), 0:max(q1))
      if(!is.null(downscaled)){
            points(q1, 
                   quantile(w, probs = seq(0.01, .99, 0.01), na.rm = T, type =4), 
                   col="blue")
      }
      par(mfrow = c(1,1)) 
}

#end    
