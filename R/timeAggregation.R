#' @title Time aggregation of the data in a field/multifield
#' @description Aggregates daily or sub-daily data annually or daily by specifiyng an aggregation function.
#' @param obj a field or multifield (can be multimember).
#' @param aggr.d a single character or a character vector indicating the daily aggregation function. 
#' Possibilities are:  "mean", "min", "max", "sum" (default is NULL).
#' @param aggr.y a single character or a character vector indicating the annual aggregation function. 
#' Possibilities are: "none", "mean", "min", "max", "sum" (default is "mean").
#' @return A field or multifield with the time dimension daily or annually aggregated (see details).
#' @details If aggr.d = NULL and the data in the field is sub-daily, annualAggregation will return an error. 
#' If the data is sub-daily and aggr.y = "none", the function will return the daily aggregated field. 
#' @author M. Iturbide \email{maibide@@gmail.com}
#' @export


timeAggregation <- function(obj, aggr.d = NULL, aggr.y = "mean"){
      
      
      
      x.coord <- getCoordinates(obj)$x
      y.coord <- getCoordinates(obj)$y
      coords <- expand.grid(1:length(x.coord), 1:length(y.coord))
      

      
      dimNames <- attr(obj$Data, "dimensions")
      if (any(grepl("var", dimNames))) {
            aux.dates <- as.POSIXlt(obj$Dates[[1]]$start)
            aux.dates.end <- as.POSIXlt(obj$Dates[[1]]$end)
            
      } else {
            aux.dates <- as.POSIXlt(obj$Dates$start)
            aux.dates.end <- as.POSIXlt(obj$Dates$end)
      }
      
      if(abs(difftime(aux.dates[1], aux.dates[2], units = "days")) < 1){
            if(is.null(aggr.d)){
                  stop("Data is sub-daily: argument aggr.d needs to be defined")
            }
            
            yrs <- aux.dates$year + 1900
            #     
            yy <- unique(yrs)
            
            dys <- lapply(1:length(yy), function(i){
                  if(i == 1){
                        aux.dates$yday[which(yrs == yy[i])] 
                  }else{
                        aux.dates$yday[which(yrs == yy[i])] + (i-1)*365
                  }
            })
            period.id <- do.call("abind", dys)
            #       if (any(dys != unique(dys)){
            
             
      
 
            if(any(attr(obj$Data, "dimensions")=="var")){
                  vars <- obj$Variable$varName
                  if(length(aggr.d)==1){
                        aggr.d <- rep(aggr.d, length(vars))    
                  }else if(length(aggr.d)!=length(vars)){
                        message("NOTE: Argument 'aggr.d' is not of the same length as the number of variables in the multifield, 
                                only the first function will be applied to all the variables")
                        aggr.d <- rep(aggr.d[1], length(vars)) 
                        
                        
                  }else{
                        aggr.d <- aggr.d   
                  }
                 
                  if (any(attr(obj$Data, "dimensions")=="member")){
                        mind <- which(attr(obj$Data, "dimensions")=="member")
                        nmem <- dim(obj$Data)[mind]
                        fin <- array(dim = c(1,nmem,length(unique(period.id)), length(x.coord), ncol = length(y.coord)))
                        lfin <- lapply(1:length(aggr.d), function(i){
                              message("Applying ", aggr.d[i], " daily aggregation function to var ", vars[i])
                              fin[,,,,] <- annualAggregationMain(obj$Data[i,,,,], FUN = aggr.d[i], id = period.id)
                              return(fin)
                        })
                        
                  }else{
                        fin <- array(dim = c(1,length(unique(period.id)), length(x.coord), ncol = length(y.coord)))  
                        lfin <- lapply(1:length(aggr.d), function(i){ 
                              message("Applying ", aggr.d[i], " daily aggregation function to var ", vars[i])
                              fin[,,,] <- annualAggregationMain(obj$Data[i,,,], FUN = aggr.d[i], id = period.id)
                              return(fin)
                        })
                  }
                  
                  
                  dato <- do.call("abind", c(lfin, along = 1))
                  
            }else{
                  message("Applying ", aggr.d, " daily aggregation function")
                  dato <- annualAggregationMain(obj$Data, FUN =aggr.d, id = period.id)
            }
            dims <- attr(obj$Data, "dimensions")
            obj$Data <- unname(dato)
            attr(obj$Data, "dimensions") <- dims
            attr(obj$Variable, "daily_agg_cellfun") <- aggr.d
            if(any(attr(obj$Data, "dimensions")=="var")){
                  obj$Dates <- list("start" = unname(tapply(obj$Dates[[1]]$start, INDEX = period.id, FUN = min)), 
                                    "end" = unname(tapply(obj$Dates[[1]]$end, INDEX = period.id, FUN = min)))       
            }else{
                  obj$Dates <- list("start" = unname(tapply(obj$Dates$start, INDEX = period.id, FUN = min)), 
                                    "end" = unname(tapply(obj$Dates$end, INDEX = period.id, FUN = min))) 
            }
      }else if(!is.null(aggr.d)){
            message("Data is not sub-daily: aggr.d will be ignored")
      }
      
      
      ####
      if(!is.null(aggr.d) & aggr.y != "none" ){
            period.id <- unname(getYearsAsINDEX(obj))
      
            if(any(attr(obj$Data, "dimensions")=="var")){
                  vars <- obj$Variable$varName
                  if(length(aggr.y)==1){
                        aggr.f <- rep(aggr.y, length(vars))    
                  }else if(length(aggr.y)!=length(vars)){
                        message("NOTE: Argument 'aggr.y' is not of the same length as the number of variables in the multifield, 
                          only the first function will be applied to all the variables")
                        aggr.f <- rep(aggr.y[1], length(vars)) 
                  
                  
                  }else{
                        aggr.f <- aggr.y   
                  }
            ###
                  if (any(attr(obj$Data, "dimensions")=="member")){
                        mind <- which(attr(obj$Data, "dimensions")=="member")
                        nmem <- dim(obj$Data)[mind]
                        fin <- array(dim = c(1,nmem,length(unique(period.id)), length(y.coord), ncol = length(x.coord)))
                        lfin <- lapply(1:length(aggr.f), function(i){
                              message("Applying ", aggr.f[i], " annual aggregation function to var ", vars[i])
                              fin[,,,,] <- annualAggregationMain(obj$Data[i,,,,], FUN = aggr.f[i], id = period.id)                                                                
                              return(fin)
                        })
                 
                  }else{
                        fin <- array(dim = c(1,length(unique(period.id)), length(y.coord), ncol = length(x.coord)))  
                        lfin <- lapply(1:length(aggr.f), function(i){ 
                              message("Applying ", aggr.f[i], " annual aggregation function to var ", vars[i])
                              fin[,,,] <- annualAggregationMain(obj$Data[i,,,], FUN = aggr.f[i], id = period.id)
                              return(fin)
                        })
                  }
                    
            
                  dato <- do.call("abind", c(lfin, along = 1))
         
            }else{
                  aggr.f <- aggr.y   
                  message("Applying ", aggr.f, " annual aggregation function")
                  dato <- annualAggregationMain(obj$Data, FUN =aggr.y, id = period.id)
            }
      
      
      
      dims <- attr(obj$Data, "dimensions")
      obj$Data <- unname(dato)
      attr(obj$Data, "dimensions") <- dims
      attr(obj$Variable, "annually_agg_cellfun") <- aggr.f
            if(any(attr(obj$Data, "dimensions")=="var")){
                  obj$Dates <- list("start" = unname(tapply(obj$Dates[[1]]$start, INDEX = period.id, FUN = min)), 
                                    "end" = unname(tapply(obj$Dates[[1]]$end, INDEX = period.id, FUN = min)))       
            }else{
                  obj$Dates <- list("start" = unname(tapply(obj$Dates$start, INDEX = period.id, FUN = min)), 
                        "end" = unname(tapply(obj$Dates$end, INDEX = period.id, FUN = min))) 
            }
      }
return(obj)
}

#end    


#' @title Annual aggregation of the data in a field
#' @description Aggregates data sub-annual data annually by specifiyng an aggregation function
#' @param data array (can be multimember).
#' @param FUN a character indicating the aggregation function. Possibilities are: "mean", "min", "max", "sum".
#' @param id as argument INDEX in function tapply. 
#' @return array  with the time dimension annually aggregated.
#' @author M. Iturbide \email{maibide@@gmail.com}

annualAggregationMain <- function(data, FUN = c("mean", "min", "max", "sum"), id){
   
      aggr <- match.arg(FUN, choices = c("none","mean", "min", "max", "sum"))


            if (any(attr(obj$Data, "dimensions")=="member")){
                  mind <- which(attr(obj$Data, "dimensions")=="member")
                  nmem <- dim(obj$Data)[mind]
                  r <- array(dim = c(1,length(unique(id)), length(y.coord), ncol = length(x.coord)))
                  
                  w <- lapply(1:nmem, function(m){
                        message("[", Sys.time(), "] Aggregating member ", m, " out of ", nmem)
                        for (i in 1:nrow(coords)){
                  
                              r[,,coords[i,2],coords[i,1]] <- tapply(data[m,, coords[i,2], coords[i,1]],
                                     INDEX = id, FUN = aggr)
                         }
                        return(r)
                  })
                  message("[", Sys.time(), "] Done.")
            dat <- unname(do.call("abind", c(w, along=mind)))
            }else{
                  dat <- array(dim = c(length(unique(id)), length(y.coord), ncol = length(x.coord)))
                  message("[", Sys.time(), "] Aggregating...")
                  for (i in 1:nrow(coords)){
                        
                        dat[,coords[i,2],coords[i,1]] <- tapply(data[, coords[i,2], coords[i,1]],
                                                             INDEX = id, FUN = aggr)
                  }
                  message("[", Sys.time(), "] Done.")
            }
      
      return(dat)
}

#end    
