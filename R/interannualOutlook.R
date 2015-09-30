#' @title Plots annual series of the inputs and the output of downscaling
#' @description Plots annual series of the inputs and the output of downscaling
#' @param obs Field of observations. 
#' @param pred Field of model data.
#' @param downscaled Field of downscaling output.  
#' @param location coordinates of a location in the geographic domain of the field.  
#' @param yrange equivalent to ylim argument for plot function  
#' 
#' @return A plot
#' @author M. Iturbide \email{maibide@@gmail.com}
#' @export

interannualOutlook <- function(obs, pred, downscaled, location = c(-40, -5), yrange = NULL){
  
  test.id2 <- getYearsAsINDEX(downscaled)
  period.id <- (getYearsAsINDEX(pred))
  train.id <- period.id[which(is.na(match(period.id, test.id2)))]
  test.id <- period.id[which(is.na(match(period.id, train.id)))]
  
  period <- unique(period.id)
  test.period <- unique(test.id)
  train.period <- period[which(is.na(match(period,test.period)))]
  
  pred.test <- subsetField(pred, years = test.period)
  pred.train <- subsetField(pred, years = train.period)
  
  obs.train <- subsetField(obs, years = train.period)
  # as reference to validate results, the subset for the 
  # observations in the test period:
  obs.test <- subsetField(obs, years = test.period) 
  
  
  
  
  # WFDEI Mean tp climatology (1981-2010) interpolated to the grid of NCEP
  #   obs.mat <- apply(obs$Data, MARGIN = 2:3, FUN = "mean")
  #   require(fields) # Used for plotting map
  x.coord <- getCoordinates(obs)$x
  y.coord <- getCoordinates(obs)$y
  #   image.plot(x.coord, y.coord, t(obs.mat), asp = 1, main = "Mean of the OBSERVATIONS in the in the full period")
  #   world(add = TRUE)
  #   points(x.coord[40], y.coord[30], cex = 2, pch = 17)
  #   
  #   # Find the closest cell in NCEP to the selected geographic point
  #   
  xo <- findInterval(location[1], x.coord)
  yo <- findInterval(location[2], y.coord)
  
  xi <- findInterval(x.coord[xo], getCoordinates(pred)$x)
  yi <- findInterval(y.coord[yo], getCoordinates(pred)$y)
  #   
  
  if (any(attr(pred$Data, "dimensions")=="member")){
    nmem <- dim(pred$Data)[which(attr(pred$Data, "dimensions")=="member")]
    #Daily time series plot
    
    ##Compute statistic
    x <- tapply(obs$Data[, yo,xo], INDEX = period.id, FUN = mean)
    
    yl <- lapply(1:nmem, function(m){
      tapply(pred$Data[m,,yi,xi], INDEX = period.id, FUN = mean)})
    y <- Reduce("+", yl)/length(yl)
    ys <- apply(simplify2array(yl), 1, sd)
    
    ytl <- lapply(1:nmem, function(m){
      tapply(pred.test$Data[m,,yi,xi], INDEX = test.id, FUN = mean)})
    yt <- Reduce("+", ytl)/length(ytl)
    
    wl <- lapply(1:nmem, function(m){
      tapply(downscaled$Data[m,,yo,xo], INDEX = test.id2, FUN = mean)})
    
    w <- Reduce("+", wl)/length(wl)
    ws <- apply(simplify2array(wl), 1, sd)
    
    
    ## ACCURACY
    ### Spearman correlation rho and the Root Mean Square Error (RMSE) as accuracy measures 
    ### for the direct and calibrated simulation in the TEST PERIOD
    
    xt <- tapply(obs.test$Data[, yo,xo], INDEX = test.id, FUN = mean)
    
    
    rmse.down <- sqrt(mean((xt - w)^2))
    bias.down <-  sum(w - xt)/sum(x)
    rho.down <- cor(x = xt, y = w, method = "spearman")
    
    rmse.direct <- sqrt(mean((xt - yt)^2))
    bias.direct <-  sum(yt - xt)/sum(x)
    rho.direct <- cor(x = xt, y = yt, method = "spearman")
    
    
    if (is.null(yrange)) {
      mi <- floor(min(c(x,y,w)))-1
      ma <-  floor(max(c(x,y,w)))
      yrange <- c(mi, ma + (ma-mi))}
    
    plot(1:length(period), x, xlim = c(0,length(period)), ylim = yrange, xlab="", xaxt = "n", 
         ylab = "Annual/seasonal mean value", cex = .6, col = NULL)
    tck <- axis(1, at = 1:length(period), labels=FALSE)
    text(tck,  par("usr")[3] - 1, xpd = TRUE, labels = (1981:2010), 
         srt = 90, cex =.6)
    #plot the sd (shadows)
    polygon(x = c(1:length(period), length(period):1), y =c (ys+y,rev(y-ys)), col = rgb(1,0,0,0.2), border = NA)
    polygon(x = c((length(train.period)+1):length(period), length(period):(length(train.period)+1)), 
            y =c (ws+w,rev(w-ws)), col = rgb(0,0,1,0.2), border = NA)
    #plot the mean (lines)
    lines(1:(length(train.period)+1), x[1:(length(train.period)+1)], lwd = 2, lty = 4, xlim = c(0,length(period)))
    lines((1+length(train.period)):length(period) , x[(1+length(train.period)):length(period)], 
          lwd = 2)
    
    lines(1:(length(train.period)+1), y[1:(length(train.period)+1)], lwd = 2, lty = 4, xlim = c(0,length(period)),
          col="red")
    lines((1+length(train.period)):length(period) , y[(1+length(train.period)):length(period)], 
          lwd = 2, col ="red")
    lines((length(train.period)+1):length(period), w, col = "blue", lwd = 2)
    
    legend(0, yrange[2] , legend = c("obs",  
                                     paste("direct: rho=", 
                                           round(rho.direct, digits = 2), ", bias=", 
                                           round(bias.direct, digits = 2), sep = ""), 
                                     paste("downscaled: rho=", 
                                           round(rho.down, digits = 2), ", bias=", 
                                           round(bias.down, digits = 2), sep = "")), 
           fill = c("black", "red", "blue"), box.lwd = 0, cex = .8)
    
    
    
    
  }else{
    #Daily time series plot
    
    
    ##Compute statistic
    x <- tapply(obs$Data[, yo,xo], INDEX = period.id, FUN = mean)
    
    y <- tapply(pred$Data[,yi,xi], INDEX = period.id, FUN = mean)
    
    w <- tapply(downscaled$Data[,yo,xo], INDEX = test.id2, FUN = mean)
    
    ## ACCURACY
    ### Spearman correlation rho and the Root Mean Square Error (RMSE) as accuracy measures 
    ### for the direct and calibrated simulation in the TEST PERIOD
    
    xt <- x[(1+length(train.period)):length(period)]
    
    yt <- y[(1+length(train.period)):length(period)]
    
    rmse.down <- sqrt(mean((xt - w)^2))
    bias.down <-  sum(w - xt)/sum(x)
    rho.down <- cor(x = xt, y = w, method = "spearman")
    
    rmse.direct <- sqrt(mean((xt - yt)^2))
    bias.direct <-  sum(yt - xt)/sum(x)
    rho.direct <- cor(x = xt, y = yt, method = "spearman")
    
    if (is.null(yrange)) {
      mi <- floor(min(c(x,y,w)))-1
      ma <-  floor(max(c(x,y,w)))
      yrange <- c(mi, ma + (ma-mi))}
    
    plot(1:length(period), x, xlim = c(0,length(period)), ylim = yrange, xlab="", xaxt = "n", 
         ylab = "Annual/seasonal mean value", cex = .6, col = NULL)
    tck <- axis(1, at = 1:length(period), labels=FALSE)
    text(tck,  par("usr")[3] - 1, xpd = TRUE, labels = (1981:2010), 
         srt = 90, cex =.6)
    #plot the mean (lines)
    lines(1:(length(train.period)+1), x[1:(length(train.period)+1)], lwd = 2, lty = 4, xlim = c(0,length(period)))
    lines((1+length(train.period)):length(period) , x[(1+length(train.period)):length(period)], 
          lwd = 2)
    
    lines(1:(length(train.period)+1), y[1:(length(train.period)+1)], lwd = 2, lty = 4, xlim = c(0,length(period)),
          col="red")
    lines((1+length(train.period)):length(period) , y[(1+length(train.period)):length(period)], 
          lwd = 2, col ="red")
    lines((length(train.period)+1):length(period), w, col = "blue", lwd = 2)
    
    legend(0, yrange[2] , legend = c("obs",  
                                     paste("direct: rho=", 
                                           round(rho.direct, digits = 2), ", bias=", 
                                           round(bias.direct, digits = 2), sep = ""), 
                                     paste("downscaled: rho=", 
                                           round(rho.down, digits = 2), ", bias=", 
                                           round(bias.down, digits = 2), sep = "")), 
           fill = c("black", "red", "blue"), box.lwd = 0, cex = .8)
    
    
  }
  
  
}

# End      
