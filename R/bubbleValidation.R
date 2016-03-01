#' @title Bubble plot 
#' @description Bubble plot for the visualization of the skill of an ensemble forecast prediction. It provides a
#'  spatially-explicit representation of the skill, resolution and reliability of a probabilistic predictive system in
#'  a single map.
#' @param mm.obj A multi-member object with predictions, either a grid or a multi-member station object as a result of
#' downscaling of a forecast using station data. See details.
#' @param obs The benchmarking observations for forecast verification. 
#' @param select.year Year within the whole verification period to display the results for.
#' @param score Logical. Whether to include or not the relative operating characteristic skill score (ROCSS). See details.
#' @param size.as.probability Logical. Whether to include the tercile probabilities (magnitude proportional to bubble radius)
#'  in the plot. See details. 
#' @importFrom scales alpha
#' @importFrom verification roc.area
#' @export
#' @details 
#' For each member, the daily predictions are averaged to obtain a single seasonal forecast. The corresponding terciles 
#' for each ensemble member are then computed for the analysis period. Thus, each particular grid point, member and season,
#' are categorized into three categories (above, between or below), according to their respective climatological 
#' terciles. Then, a probabilistic forecast is computed year by year by considering the number of members falling 
#' within each category. For instance, probabilities below 1/3 are very low, indicating that a minority of the members 
#' falls in the tercile. Conversely, probabilities above 2/3 indicate a high level of member agreement (more than 66\% of members
#' falling in the same tercile). Color represents the tercile with the highest probability for the selected year. The bubble size 
#' indicates the probability of that tercile. This option is not plotted if the size.as.probability argument is FALSE.
#' 
#' Finally, the ROC Skill Score (ROCSS) is computed. For each tercile, it provides a quantitative measure of the forecast skill,
#' and it is commonly used to evaluate the performance of probabilistic systems (Joliffe and Stephenson 2003). The value of 
#' this score ranges from 1 (perfect forecast system) to -1 (perfectly bad forecast system). A value zero indicates no skill 
#' compared with a random prediction. The transparency of the bubble is associated to the ROCSS (negative values are
#' plotted with x).  This option is not plotted if the score argument is FALSE.
#' @note The computation of climatological terciles requires a representative period to obtain meaningful results.
#' @author J. Bedia, M.D. Frias and J. Fernandez 
#' @family visualization
#' @references
#'  Jolliffe, I. T. and Stephenson, D. B. 2003. Forecast Verification: A Practitioner's Guide in 
#'  Atmospheric Science, Wiley, NY


bubbleValidation <- function(mm.obj, obs, select.year, score = TRUE, size.as.probability = TRUE) {
      # Transform both data to the same grid (model grid)
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
      if (!identical(c("time", "lat", "lon"), obs.dimNames)) {
            stop("The observed reference must be a 3D array of the form [time,lat,lon]")
      }
      obs <- interpData(obs, new.Coordinates = getGrid(mm.obj), method = "nearest")  
      x.mm <- mm.obj$xyCoords$x
      y.mm <- mm.obj$xyCoords$y
      yrs <- getYearsAsINDEX(mm.obj)
      yy <- unique(yrs)
      if (!select.year %in% yy) {
            stop("Target year outside temporal data range")
      }
      iyear <- which(yy[1]:yy[length(yy)] == select.year)
      lon.dim <- grep("lon", mm.dimNames)
      lat.dim <- grep("lat", mm.dimNames)
      member.dim <- grep("member", mm.dimNames)
      n.mem <- dim(mm.obj$Data)[member.dim]
      # Computation of terciles and exceedance probabilities
      # yearmean changes the data dimension. time.dim is in the first dimension!!
      y.mean <- apply(mm.obj$Data, MARGIN = c(lat.dim, lon.dim, member.dim), FUN = function(x) {
            tapply(x, INDEX = yrs, FUN = mean, na.rm = TRUE)
      })
      terciles <- apply(y.mean, MARGIN = c(2, 3, 4), FUN = quantile, c(1/3, 2/3), na.rm = TRUE)
      # Compute the probability for each tercile
      t.u <- array(dim = dim(y.mean)[1:3])
      t.l <- t.u
      t.m <- t.u
      prob <- array(dim = c(3,dim(y.mean)[1:3]))
      for (i in seq.int(1, dim(y.mean)[1])) {
            t.u[i, , ] <- apply(y.mean[i, , , ] > terciles[2, , , ], MARGIN = c(1, 2), FUN = sum, na.rm = TRUE) / n.mem
            t.l[i, , ] <- apply(y.mean[i, , , ] < terciles[1, , , ], MARGIN = c(1, 2), FUN = sum, na.rm = TRUE) / n.mem
            t.m[i, , ] <- 1 - t.u[i, , ] - t.l[i, , ]
            prob[1, i, , ] <- t.l[i, , ]
            prob[2, i, , ] <- t.m[i, , ]
            prob[3, i, , ] <- t.u[i, , ]
      }
      # Maximum probability from the terciles
      max.prob <- apply(prob, MARGIN = c(2, 3, 4), FUN = max, na.rm = TRUE)
      # Tercile for the maximum probability from the terciles
      t.max.prob <- apply(prob, MARGIN = c(2, 3, 4), FUN = which.max)
      # Terciles for the observations
      obs.y.mean <- apply(obs$Data, MARGIN = c(2, 3), FUN = function(x) {
            tapply(x, INDEX = yrs, FUN = mean, na.rm = TRUE)
      })
      obs.terciles <- apply(obs.y.mean, MARGIN = c(2, 3), FUN = quantile, c(1/3, 2/3), na.rm = TRUE)    
      obs.t.u <- array(dim = dim(obs.y.mean))
      obs.t.l <- obs.t.u
      obs.t.m <- obs.t.u
      for (i in seq.int(1, dim(obs.y.mean)[1])) {
            obs.t.u[i, , ] <- (obs.y.mean[i, , ] > obs.terciles[2, , ])
            obs.t.l[i, , ] <- (obs.y.mean[i, , ] < obs.terciles[1, , ])
            obs.t.m[i, , ] <- (obs.y.mean[i, , ] >= obs.terciles[1, , ] & obs.y.mean[i, , ] <= obs.terciles[2, , ])
      }
      obs.t <- obs.t.u * 1 + obs.t.l * -1 # 1 upper tercile, 0 middle tercile, -1 lower tercile
      # Filter points with observations in model data 
      # Select a year and eliminate de NaN cases detected for the observations. 
      v.max.prob <- as.vector(max.prob[iyear, , ])
      v.t.max.prob <- as.vector(t.max.prob[iyear, , ])
      v.nans <- complete.cases(as.vector(obs.t[iyear, , ]))
      ve.max.prob <- v.max.prob[v.nans]
      if (!size.as.probability) {
            ve.max.prob <- rep(0.5, length(ve.max.prob))
      }
      df <- data.frame(max.prob = ve.max.prob, t.max.prob = v.t.max.prob[v.nans])
      df$color <- "black"
      df$color[df$t.max.prob == 3] <- "red"
      df$color[df$t.max.prob == 2] <- "gold"
      df$color[df$t.max.prob == 1] <- "blue"
      yx <- as.matrix(expand.grid(y.mm, x.mm))
      nn.yx <- yx[v.nans, ]
      score = TRUE
      if (score) { # Compute ROCSS
            v.score <- c()
            count <- 1
            for (ilon in seq(1,length(x.mm))) {
                  for (ilat in seq(1,length(y.mm))) {
                        n.nan <- sum(is.na(obs.t[,ilat,ilon]))
                        if (n.nan == 0) {
                              # Compute the score only for the tercile with maximum probability
                              select.tercile <- t.max.prob[iyear,ilat,ilon] 
                              if (select.tercile == 1) {
                                    res <- suppressWarnings(roc.area(obs.t.l[ , ilat, ilon], t.l[ , ilat, ilon]))
                                    v.score[count] <- res$A*2 - 1
                              } 
                              if (select.tercile == 2) {
                                    res <- suppressWarnings(roc.area(obs.t.m[ ,ilat, ilon], t.m[ , ilat, ilon]))
                                    v.score[count] <- res$A*2 - 1
                              }
                              if (select.tercile == 3) {
                                    res <- suppressWarnings(roc.area(obs.t.u[,ilat,ilon], t.u[,ilat,ilon]))
                                    v.score[count] <- res$A*2 - 1
                              }      
                              count <- count + 1
                        }
                  }  
            }
            # Select positive score values from negative values
            pos.val <- v.score >= 0
            neg.val <- v.score < 0
      }
      # Bubble plot
      par(bg = "white", mar = c(3, 3, 1, 5))
      if (score) {
            plot(nn.yx[pos.val, 2], nn.yx[pos.val, 1], cex = df$max.prob[pos.val] * 3, col = alpha(df$color[pos.val], v.score[pos.val]), pch = 19, xlab = "", ylab = "")
            points(nn.yx[neg.val, 2], nn.yx[neg.val, 1], pch = 4, cex = 0.75)
      } else {
            plot(nn.yx[ , 2], nn.yx[ , 1], cex = df$max.prob * 3, col = df$color, pch = 19, xlab = "", ylab = "")
      }
      draw.world.lines()
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
      if (score) {
            legend('right', c("T1", "T2", "T3", "Negat.\nScore"), pch = c(19, 19, 19, 4), col = c("blue", "gold", "red", "black"), inset = c(0,0), xpd = TRUE, bty = "n")
      } else {
            legend('right', c("T1", "T2", "T3"), pch = rep(19, 3), col = c("blue", "gold", "red"), inset = c(0, 0), xpd = TRUE, bty = "n")
      }
}
# End

