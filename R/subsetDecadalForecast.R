#' @title Subset decadal forecast by initialization years
#' 
#' @description Subset decadal forecast by initialization years
#' 
#' @param data field obtained from the output of function loadDecadalForecast. 
#' @param ly Vector of integers indicating the lead years of initialization. 
#' Default is ly = c(1,2) (the two years before).
#' @return A field object with an additional dimension corresponding to the selected lead years. 
#' @author M. Iturbide
#' @importFrom abind abind
#' @export
#' @examples \dontrun{
#' latLim <- c(35,45)
#' lonLim <-  c(-10,2)
#' season <- 3:5
#' period <- 1981:1982
#' loginUDG(username = "myuser", password = "mypassword") #type help(loginUDG)
#' tas_decadalForecast <- loadDecadalForecast(
#'    dataset = "http://www.meteo.unican.es/tds5/dodsC/specs/gfdl_specs_decadal.ncml", 
#'    latLim = latLim, 
#'    lonLim = lonLim,
#'    var = "tas", 
#'    dictionary = F, 
#'    years = period, 
#'    season = season)
#' 
#' data(tas_decadalForecast)
#' # data corresponding to the initialitation of the two years before
#' subsetDecadalForecast(tasDECA, ly = c(1,2))
#' }

subsetDecadalForecast <- function(data, ly = c(1,2)){
      init <- data$InitializationDates      
      dates <- data$Dates$start
      m <- months(as.Date(dates))
      if (any(m == "diciembre")) {
            elim <- TRUE
            init <- init[-1]
            if (any(ly == 9)) {
            message("lead year 9 does not contain data for december and will be ignored")
                  if (length(ly) == 1) {
                        stop("please, request another lead year")
                  } else {
                        ly <- ly[-which(ly == 9)]
                  }
            }
      }
      y <- getYearsAsINDEX(data)
      yu <- unique(y)
      months <- (which(y[1:12] != y[1])[1] - 1)/2
      if (elim == T) {months <- round(months)}
      if (is.na(months) & length(yu) > 1) {
            months <- 12
      } else if (is.na(months) & length(yu) == 1) {
            months <- length(y)/(length(data$InitializationDates) - length(unique(y)) + 1)
      }
      InitializationDates <- list()
      start <- list()
      dat <- list()
      for (n in 1:length(ly)) {
            yl <- (length(data$InitializationDates) - length(unique(y)) + 1) - ly[n]
            if (elim == T) {yl <- yl - 1}
            IND <- list()
            yinit.sub <- list()
                  for (i in 1:length(yu)) {
                        year <- which(y == yu[i])
                        if (elim == T) {year <- year[-(1:(months-1))]}
                        a <- yl*months
                        ind <- lapply(1:length(a), function(x){
                        (a[x] - months + 1):a[x]
                  
                        })
            
                        yinit <- as.numeric(substr(init, 1, 4))
                        initdates <- yu[i]-(length(yinit)-length(unique(y))+1)+yl
                        initdates2 <- numeric()
                              for (k in 1:length(initdates)){
                                    initdates2[k] <- which(yinit == initdates[k])
                              }
                        yinit.sub[[i]] <- initdates2
                  
                        IND[[i]] <- year[unlist(ind)]
                         
                  }
            
        
            
            ind <- sort(unlist(IND))
            
            
            if(any(attr(data$Data, "dimensions") == "member")){
                  subs <- data$Data[,ind,,]
                  arr <- array(data = NA, dim = c(1, dim(subs)))
                  arr[1,,,,]  <- subs
                  dat[[n]] <- arr
            }else{
                  subs<- data$Data[ind,,]
                  arr <- array(data = NA, dim = c(1, dim(subs)))
                  arr[1,,,]  <- subs
                  dat[[n]] <- arr
            }
            ind.init <- sort(unique(unlist(yinit.sub)))
            if (length(ind.init ) == 0) {stop("ly (lead years) outside the extent of the number of initializations")}
            InitializationDates[[n]] <- init[ind.init]
            start <- dates[ind]
      
    
      
      }
      if(length(ly)==1){
            datbinded <- dat[[1]]
      }else{
      datbinded <- unname(abind(dat, along = 1))
      }
      names(InitializationDates) <- paste("Lead year = ", ly)
      data$InitializationDates <- InitializationDates
      data$Dates$start <- start
      data$Dates$end <- start
      
     
      attr(datbinded, "dimensions") <- c("runtime", attr(data$Data, "dimensions"))
      data$Data <- datbinded
      
      return(data)
}

#End



