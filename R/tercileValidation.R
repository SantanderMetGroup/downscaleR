#' @title Tercile plot for visualization of the skill of an ensemble forecast prediction
#' 
#' @description Tercile plot for the visualization of the skill of an ensemble forecast prediction.
#' 
#' @param mm.obj A multi-member object, either a field or a multi-member station object as a result of
#' downscaling of a forecast using station data. See details.
#' @param obs The benchmarking observations for forecast verification
#' @param stationId In case of multimember multistation objects, one station can be selected to plot
#'  the diagram. Otherwise ignored.
#' @param color.pal Color palette for the representation of the probabilities. Default to \code{"ypb"} (yellow-pink-blue),
#'  suitable for the visualization in black and white printing devices. \code{"reds"} for a white-red transition.
#' 
#' @importFrom abind asub
#' @importFrom verification roc.area
#' 
#' @export
#' 
#' @details 
#'  
#' For each member, the daily predictions are averaged to obtain a single seasonal forecast. For
#' rectangular spatial domains (i.e., for fields), the spatial average is first computed (with a warning) to obtain a
#' unique series for the whole domain. The corresponding terciles for each ensemble member are then computed
#' for the analysis period. Thus, data is converted converted to a series of tercile categories by considering values
#'  above, between or below the terciles of the whole period. The probability of a member to fall into the observed tercile
#'  is represented by the colorbar. For instance, probabilities below 1/3 are very low, indicating that a minority of the members 
#'  falls in the tercile. Conversely, probabilities above 2/3 indicate a high level of member agreement (more than 66\% of members
#'  falling in the same tercile). The observed terciles (the events that actually occurred) are represented by the white circles.
#'  
#'  Finally, the ROC Skill Score (ROCSS) is indicated in the secondary (right) Y axis. For each tercile, it provides a 
#'  quantitative measure of the forecast skill, and it is commonly used to evaluate the performance of probabilistic systems
#'  (Joliffe and Stephenson 2003). The value of this score ranges from 1 (perfect forecast system) to -1 
#'  (perfectly bad forecast system). A value zero indicates no skill compared with a random prediction.
#'  
#'  In case of multimember fields, the field is spatially averaged to obtain one single time series
#'  for each member prior to data analysis, with a warning. In case of multimember stations, one single station
#'  can be selected through the \code{stationId} argument, otherwise all station series are also averaged.
#'   
#' 
#' @note The computation of climatological terciles requires a representative period to obtain meaningful results.
#' 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}, M.D. Frias and J, Fernandez based on the original diagram 
#' conceived by A. Cofino.
#' 
#' @family visualization
#' 
#' @references
#' Diez, E., Orfila, B., Frias, M.D., Fernandez, J., Cofino, A.S., Gutierrez, J.M., 2011. 
#' Downscaling ECMWF seasonal precipitation forecasts in Europe using the RCA model.
#'  Tellus A 63, 757-762. doi:10.1111/j.1600-0870.2011.00523.x
#'    
#'  Jolliffe, I. T. and Stephenson, D. B. 2003. Forecast Verification: A Practitioner's Guide in 
#'  Atmospheric Science, Wiley, NY
#'  

tercileValidation <- function(mm.obj, obs, stationId = NULL, color.pal = c("ypb", "reds")) {
      color.pal <- match.arg(color.pal, c("ypb", "reds"))
      mm.dimNames <- attr(mm.obj$Data, "dimensions")
      obs.dimNames <- attr(obs$Data, "dimensions")
      if (!("member" %in% mm.dimNames)) {
            stop("The input data for 'multimember' is not a multimember field")
      }
      if ("member" %in% obs.dimNames) {
            stop("The verifying observations can't be a multimember")
      }
      if ("var" %in% mm.dimNames | "var" %in% obs.dimNames) {
            stop("Multifields are not a valid input")
      }
      # Preparation of a 2D array with "member" and "time" dimensions for the target
      is.mm.station <- ifelse(exists("Metadata", where = mm.obj), TRUE, FALSE)
      is.obs.station <- ifelse(exists("Metadata", where = obs), TRUE, FALSE)
      if (identical(mm.dimNames, c("member", "time", "lat", "lon"))) {
            warning("The results presented are the spatial mean of the input field")
            lat.dim.index <- grep("lat", mm.dimNames)
            lon.dim.index <- grep("lon", mm.dimNames)
            mar <- setdiff(1:length(mm.dimNames), c(lat.dim.index, lon.dim.index))
            arr <- apply(mm.obj$Data, mar, mean, na.rm = TRUE)
            attr(arr, "dimensions") <- mm.dimNames[mar]
            x.mm <- mm.obj$xyCoords$x
            y.mm <- mm.obj$xyCoords$y
      } else {
            if (identical(mm.dimNames, c("member", "time"))) {
                  arr <- mm.obj$Data
                  if (is.mm.station) {
                        x.mm <- mm.obj$xyCoords[ ,1]
                        y.mm <- mm.obj$xyCoords[ ,2]
                  } else {
                        x.mm <- mm.obj$xyCoords$x
                        y.mm <- mm.obj$xyCoords$y
                  }
            } else {
                  if (identical(mm.dimNames, c("member", "time", "station"))) {
                        if (is.null(stationId)) {
                              warning("The results presented are the mean of all input stations")
                              mar <- match(setdiff(mm.dimNames, "station"), mm.dimNames)
                              arr <- apply(mm.obj$Data, mar, mean, na.rm = TRUE)
                              attr(arr, "dimensions") <- mm.dimNames[mar]
                              x.mm <- mm.obj$xyCoords[ ,1]
                              y.mm <- mm.obj$xyCoords[ ,2]
                        } else {
                              idx <- grep(stationId, mm.obj$Metadata$station_id)
                              if (identical(idx, integer(0))) {
                                    stop("The 'stationId' provided was not found")
                              }
                              st.dim.index <- grep("station", mm.dimNames)
                              arr <- asub(mm.obj$Data, idx, st.dim.index)
                              attr(arr, "dimensions") <- setdiff(mm.dimNames, "station")
                              x.mm <- mm.obj$xyCoords[idx,1]
                              y.mm <- mm.obj$xyCoords[idx,2]
                        }
                  } else {
                        stop("Invalid input data array")
                  }
            }
      }
      # Preparation of a "time" 1D vector vec for the target 
      if (identical(obs.dimNames, c("time", "lat", "lon"))) {
            # Spatial consistency check
            x.obs <- obs$xyCoords$x
            y.obs <- obs$xyCoords$y
            lat.dim.index <- grep("lat", obs.dimNames)
            lon.dim.index <- grep("lon", obs.dimNames)
            mar <- setdiff(1:length(obs.dimNames), c(lat.dim.index, lon.dim.index))
            arr.obs <- apply(obs$Data, mar, mean, na.rm = TRUE)
            attr(arr.obs, "dimensions") <- obs.dimNames[mar]
      } else {
            if (identical(obs.dimNames, "time")) {
                  arr.obs <- obs$Data
            } else {
                  if (identical(obs.dimNames, c("time", "station"))) {
                        if (is.null(stationId)) {
                              mar <- match(setdiff(obs.dimNames, "station"), obs.dimNames)
                              arr.obs <- apply(obs$Data, mar, mean, na.rm = TRUE)
                              attr(arr.obs, "dimensions") <- obs.dimNames[mar]
                        } else {
                              idx <- grep(stationId, obs$Metadata$station_id)
                              if (identical(idx, integer(0))) {
                                    stop("The 'stationId' provided was not found in the 'obs' dataset")
                              }
                              st.dim.index <- grep("station", obs.dimNames)
                              arr.obs <- asub(obs$Data, idx, st.dim.index)
                              attr(arr.obs, "dimensions") <- setdiff(obs.dimNames, "station")
                        }
                  }
            }
      }
      # Temporal matching check (obs-pred)
      obs.dates <- as.POSIXlt(obs$Dates$start)
      mm.dates <- as.POSIXlt(mm.obj$Dates$start)
      if (!identical(obs.dates$yday, mm.dates$yday) || !identical(obs.dates$year, mm.dates$year)) {
            stop("Forecast and verifying observations are not coincident in time")
      }
      mm.dates <- NULL
      yrs <- getYearsAsINDEX(obs)
      yy <- unique(yrs)
      # Computation of terciles and exceedance probabilities
      n.mem <- dim(arr)[1]
      l <- lapply(1:n.mem, function(x) tapply(arr[x, ], INDEX = yrs, FUN = mean, na.rm = TRUE))
      aux <- do.call("rbind", l)
      terciles <- apply(aux, 1, quantile, probs = c(1/3, 2/3), na.rm = TRUE)    
      t.u <- apply(aux > terciles[2, ], 2, sum) / n.mem
      t.l <- apply(aux < terciles[1, ], 2, sum) / n.mem
      t.m <- 1-t.u-t.l
      cofinogram.data <- cbind(t.l,t.m,t.u)
      # Benchmark
      obs.mean <- tapply(arr.obs, yrs, mean, na.rm = TRUE)
      obs.terciles <- quantile(obs.mean, probs = c(1/3, 2/3), na.rm = TRUE)
      obs.t.u <- obs.mean > obs.terciles[2]
      obs.t.l <- obs.mean < obs.terciles[1]
      obs.t.m <- obs.mean >= obs.terciles[1] & obs.mean <= obs.terciles[2]
      obs.t <- obs.t.u*1+obs.t.l*-1
      # Color selection
      cbar <- switch(color.pal,
               "ypb" = rev(c("#000066FF","#0000C1FF","#1600FFFF","#5D00FFFF","#A412EDFF","#EB3FC0FF","#FF6D92FF","#FF9A65FF","#FFC738FF","#FFF50AFF")),
               "reds" = c("#FFFFFF","#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C", "#CB181D","#A50F15","#67000D"))
      brks <- c(seq(0,1,length=length(cbar)+1))
      par(oma = c(0, 0, 0, 6))
      image(yy, c(-1,0,1), cofinogram.data, breaks=brks, col=cbar, ylab="", xlab="", asp = 1, yaxt="n", bty = "n", axes = FALSE)
      axis(1, at = yy)      
      points(yy, obs.t, pch = 21, bg = "white")
      axis(2, at=-1:1, labels=c("Below", "Normal", "Above"), las="2")
      # Area underneath a ROC curve
      roca.t.u <- suppressWarnings(roc.area(obs.t.u, t.u))
      roca.t.l <- suppressWarnings(roc.area(obs.t.l, t.l))
      roca.t.m <- suppressWarnings(roc.area(obs.t.m, t.m))
      # ROCSS
      rocss.t.u <- roca.t.u$A*2-1
      rocss.t.l <- roca.t.l$A*2-1
      rocss.t.m <- roca.t.m$A*2-1
      # Add skill score values to the plot
      axis(4, at=-1:1, labels=c(round(rocss.t.l,2), round(rocss.t.m,2), round(rocss.t.u,2)), las="2")
      par(oma = c(0, 0, 0, 3))
      image.plot(add = TRUE, legend.only = TRUE, breaks = brks, col = cbar, smallplot = c(0.96,0.99,0.2,0.8), zlim=c(0,1), legend.lab="Prob")
      mtext("ROCSS", side=3, line=-1.5, adj=1, font=2)
}
# End
