# Internal function to compute the time bounds of each verification time
timeBounds <- function(dic, foreDates) { 
      if (!is.null(dic)) {
            if (isTRUE(dic$doDailyMean)) {
                  foreDates <- foreDates[seq.int(1, length(foreDates), 4)]
                  ltb <- as.difftime(0, format = "%H", units = "hours")
                  utb <- as.difftime(24, format = "%H", units = "hours")
            } else {
                  ltb <- as.difftime(dic$lower_time_bound, format = "%H", units = "hours")
                  utb <- as.difftime(dic$upper_time_bound, format = "%H", units = "hours")
            }
            foreDatesStart <- as.POSIXlt(foreDates + ltb)
            foreDatesEnd <- as.POSIXlt(foreDates + utb)
            foreDates <- foreDatesStart
      } else {
            varTimeStep <- difftime(foreDates[2], foreDates[1])
            foreDatesEnd <- foreDates + varTimeStep
      }
      foreDatesEnd <- as.POSIXlt(foreDatesEnd)
      return(list("start" = foreDates, "end" = foreDatesEnd))
}
# End