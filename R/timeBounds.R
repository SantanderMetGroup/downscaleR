# Internal function to compute the time bounds of each verification time
timeBounds <- function(dic, foreDates) { 
      foreDates <- as.POSIXlt(foreDates, tz = "GMT")
      if (!is.null(dic)) {
            if (!is.na(dic$dailyAggr)) {
                  foreDates <- foreDates[seq.int(1, length(foreDates), 4)]
                  ltb <- as.difftime(0, format = "%H", units = "hours")
                  utb <- as.difftime(24, format = "%H", units = "hours")
            } else {
                  ltb <- as.difftime(dic$lower_time_bound, format = "%H", units = "hours")
                  utb <- as.difftime(dic$upper_time_bound, format = "%H", units = "hours")
            }
            foreDatesStart <- as.POSIXlt(foreDates + ltb, tz = "GMT")
            foreDatesEnd <- as.POSIXlt(foreDates + utb, tz = "GMT")
            foreDates <- foreDatesStart
      } else {
            varTimeStep <- difftime(foreDates[2], foreDates[1])
            foreDatesEnd <- foreDates + varTimeStep
      }
      return(list("start" = format(as.POSIXct(foreDates, tz = "GMT"), format = "%Y-%m-%d %H:%M:%S", usetz = TRUE), "end" = format(as.POSIXct(foreDatesEnd, tz = "GMT"), format = "%Y-%m-%d %H:%M:%S", usetz = TRUE)))
}
# End


