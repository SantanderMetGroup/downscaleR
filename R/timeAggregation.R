#' @title Time aggregation of the data in a field/multifield
#' @description Aggregates data daily, monthly or annually by specifiyng an aggregation function.
#' @param obj a field or multifield (can be multimember).
#' @param aggr.d a single character or a character string indicating the daily aggregation function. 
#' Possibilities are:  "mean", "min", "max", "sum" (default is NULL).
#' @param aggr.m a single character or a character string indicating the monthly aggregation function. 
#' Possibilities are:  "mean", "min", "max", "sum" (default is NULL).
#' @param aggr.y a single character or a character string indicating the annual aggregation function. 
#' Possibilities are: "mean", "min", "max", "sum" (default is "NULL").
#' @return A field or multifield with the time dimension daily, monthly or annually aggregated (see details).
#' @details To annually or monthly aggregate data, aggr.d and/or aggr.m might be also specified to apply different aggregations 
#' for each time scale. If the data in the field is sub-daily and aggr.d is not specified, a WARNING message will be returned.
#' @author M. Iturbide \email{maibide@@gmail.com}
#' @export



timeAggregation <- function(obj, aggr.d = NULL, aggr.m = NULL, aggr.y = NULL){
      
      
      
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
      
      #------------------daily-----------------------------------
      
      if(abs(difftime(aux.dates[1], aux.dates[2], units = "days")) < 1){
            if(is.null(aggr.d)){
                  warning("Data is sub-daily, define aggr.d for daily aggregation")
            }else{
            
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
                              fin[,,,,] <- annualAggregationMain(obj$Data[i,,,,], FUN = aggr.d[i], id = period.id, y.coord, x.coord, coords, obj)
                              return(fin)
                        })
                        
                  }else{
                        fin <- array(dim = c(1,length(unique(period.id)), length(x.coord), ncol = length(y.coord)))  
                        lfin <- lapply(1:length(aggr.d), function(i){ 
                              message("Applying ", aggr.d[i], " daily aggregation function to var ", vars[i])
                              fin[,,,] <- annualAggregationMain(obj$Data[i,,,], FUN = aggr.d[i], id = period.id, y.coord, x.coord, coords, obj)
                              return(fin)
                        })
                  }
                  
                  
                  dato <- do.call("abind", c(lfin, along = 1))
                  
            }else{
                  message("Applying ", aggr.d, " daily aggregation function")
                  dato <- annualAggregationMain(obj$Data, FUN =aggr.d, id = period.id, y.coord, x.coord, coords, obj)
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
            }
            
      }else if(!is.null(aggr.d)){
            message("Data is not sub-daily: aggr.d will be ignored")
      }
      
      ###
      
      if(!is.null(aggr.m)){
            
            if(abs(difftime(aux.dates[1], aux.dates[2], units = "days")) > 7){
                  message("Data is not sub-monthly: aggr.m will be ignored")
            }else{
                  
                  if (any(grepl("var", dimNames))) {
                        aux.dates <- as.POSIXlt(obj$Dates[[1]]$start)
                        aux.dates.end <- as.POSIXlt(obj$Dates[[1]]$end)
                        
                  } else {
                        aux.dates <- as.POSIXlt(obj$Dates$start)
                        aux.dates.end <- as.POSIXlt(obj$Dates$end)
                  }
                  
                  yrs <- aux.dates$year + 1900
                  #     
                  yy <- unique(yrs)
                  
                  mys <- lapply(1:length(yy), function(i){
                        if(i == 1){
                              aux.dates$mon[which(yrs == yy[i])] 
                        }else{
                              aux.dates$mon[which(yrs == yy[i])] + (i-1)*12
                        }
                  })
                  period.id <- do.call("abind", mys)
                  #       if (any(dys != unique(dys)){
                  
                  
                  
                  #---------------monthly--------------------------------------
                  if(any(attr(obj$Data, "dimensions")=="var")){
                        vars <- obj$Variable$varName
                        if(length(aggr.m)==1){
                              aggr.m <- rep(aggr.m, length(vars))    
                        }else if(length(aggr.m)!=length(vars)){
                              message("NOTE: Argument 'aggr.m' is not of the same length as the number of variables in the multifield, 
                                      only the first function will be applied to all the variables")
                              aggr.m <- rep(aggr.m[1], length(vars)) 
                              
                              
                        }else{
                              aggr.m <- aggr.m   
                        }
                        
                        if (any(attr(obj$Data, "dimensions")=="member")){
                              mind <- which(attr(obj$Data, "dimensions")=="member")
                              nmem <- dim(obj$Data)[mind]
                              fin <- array(dim = c(1,nmem,length(unique(period.id)), length(x.coord), ncol = length(y.coord)))
                              lfin <- lapply(1:length(aggr.m), function(i){
                                    message("Applying ", aggr.m[i], " monthly aggregation function to var ", vars[i])
                                    fin[,,,,] <- annualAggregationMain(obj$Data[i,,,,], FUN = aggr.m[i], id = period.id, y.coord, x.coord, coords, obj)
                                    return(fin)
                              })
                              
                        }else{
                              fin <- array(dim = c(1,length(unique(period.id)), length(x.coord), ncol = length(y.coord)))  
                              lfin <- lapply(1:length(aggr.m), function(i){ 
                                    message("Applying ", aggr.m[i], " monthly aggregation function to var ", vars[i])
                                    fin[,,,] <- annualAggregationMain(obj$Data[i,,,], FUN = aggr.m[i], id = period.id, y.coord, x.coord, coords, obj)
                                    return(fin)
                              })
                        }
                        
                        
                        dato <- do.call("abind", c(lfin, along = 1))
                        
                  }else{
                        message("Applying ", aggr.m, " monthly aggregation function")
                        dato <- annualAggregationMain(obj$Data, FUN =aggr.m, id = period.id, y.coord, x.coord, coords, obj)
                  }
                  dims <- attr(obj$Data, "dimensions")
                  obj$Data <- unname(dato)
                  attr(obj$Data, "dimensions") <- dims
                  attr(obj$Variable, "monthly_agg_cellfun") <- aggr.m
                  if(any(attr(obj$Data, "dimensions")=="var")){
                        obj$Dates <- list("start" = unname(tapply(obj$Dates[[1]]$start, INDEX = period.id, FUN = min)), 
                                          "end" = unname(tapply(obj$Dates[[1]]$end, INDEX = period.id, FUN = min)))       
                  }else{
                        obj$Dates <- list("start" = unname(tapply(obj$Dates$start, INDEX = period.id, FUN = min)), 
                                          "end" = unname(tapply(obj$Dates$end, INDEX = period.id, FUN = min))) 
                  }
                  
            }
      }
      
      #----------------------annually-------------------------------
      
      ####
      if(!is.null(aggr.y)){
            
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
                              fin[,,,,] <- annualAggregationMain(obj$Data[i,,,,], FUN = aggr.f[i], id = period.id, y.coord, x.coord, coords, obj)                                                                
                              return(fin)
                        })
                        
                  }else{
                        fin <- array(dim = c(1,length(unique(period.id)), length(y.coord), ncol = length(x.coord)))  
                        lfin <- lapply(1:length(aggr.f), function(i){ 
                              message("Applying ", aggr.f[i], " annual aggregation function to var ", vars[i])
                              fin[,,,] <- annualAggregationMain(obj$Data[i,,,], FUN = aggr.f[i], id = period.id, y.coord, x.coord, coords, obj)
                              return(fin)
                        })
                  }
                  
                  
                  dato <- do.call("abind", c(lfin, along = 1))
                  
            }else{
                  aggr.f <- aggr.y   
                  message("Applying ", aggr.f, " annual aggregation function")
                  dato <- annualAggregationMain(obj$Data, FUN =aggr.y, id = period.id, y.coord, x.coord, coords, obj)
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
            
      }else if(is.null(c(aggr.d, aggr.m, aggr.y))){
            stop("Undefined aggregation function")
      } 
      return(obj)
}

#end    


#' @title Time aggregation of the data in a field
#' @description Aggregates data daily, monthly or annually by specifiyng an aggregation function
#' @param data array (can be multimember).
#' @param FUN a character indicating the aggregation function. Possibilities are: "mean", "min", "max", "sum".
#' @param id as argument INDEX in function tapply. 
#' @param y.coord y.coord
#' @param x.coord x.coord
#' @param coords coords
#' @param obj obj 
#' @return array  with the time dimension aggregated.
#' @author M. Iturbide \email{maibide@@gmail.com}

annualAggregationMain <- function(data, FUN = c("mean", "min", "max", "sum"), id, y.coord, x.coord, coords, obj){
   
      aggr <- match.arg(FUN, choices = c("mean", "min", "max", "sum"))


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
