#' @title Extract field data from PCs and EOFs
#' @description Extracts field data from PCs and EOFs
#' @param pc.obj a list returned by function \code{\link{prinComp}}. Default is "NULL" (see Details).
#' @param orig logical (default is FALSE). If TRUE and argument 'pc.obj' is used, 
#' field data is reovered from the 'orig' slot, that contains the original 
#' field data before PCs calculation with function \code{\link{prinComp}}. 
#' @param EOFs matrix with EOFs in rows (default is NULL).
#' @param PCs  matrix with PCs in columns (default is "NULL").
#' @details The field/multifield can be built from the EOFs and PCs in 'pc.obj' or from matrixes of EOFs and PCs.   
#' @return A field/multifield (if 'pc.obj' is used) or a 2D matrix. 
#' @author M. Iturbide \email{maibide@@gmail.com}

PC2field <- function(pc.obj = NULL, orig = FALSE, EOFs = NULL, PCs = NULL){
      if(!is.null(pc.obj)){
            nvars <- length(pc.obj)-1
            pc.obj$Variable$varName <- names(pc.obj)[-length(pc.obj)]
            
#members remaining             
            if (length(pc.obj)>2){ 
                  pc.obj$Data <- array(data = NA, dim = c(nvars, length(attr(pc.obj, "dates_start")), 
                                                          length(attr(pc.obj, "yCoord")), length(attr(pc.obj, "xCoord"))))
                  for (i in 1:nvars){
                        if(orig == FALSE){
                              pc.obj$Data[i,,,] <- mat2Dto3Darray(
                                    t(pc.obj[[i]][[1]]$EOFs %*% t(pc.obj[[i]][[1]]$PCs)),
                                    attr(pc.obj, "xCoord"), attr(pc.obj, "yCoord")
                              )
                        }else{
                                    
                              pc.obj$Data[i,,,] <- mat2Dto3Darray(
                                    pc.obj[[i]][[1]]$orig,
                                    attr(pc.obj, "xCoord"), attr(pc.obj, "yCoord")
                              )
                        }
                  }
                  attr(pc.obj$Data, "dimensions") = c("var", "time", "lat", "lon")
            }else{
                  pc.obj$Data <- array(data = NA, dim = c(length(attr(pc.obj, "dates_start")), 
                                                          length(attr(pc.obj, "yCoord")), length(attr(pc.obj, "xCoord"))))

                        if(orig == FALSE){
                              pc.obj$Data[,,] <- mat2Dto3Darray(
                                    t(pc.obj[[1]][[1]]$EOFs %*% t(pc.obj[[1]][[1]]$PCs)),
                                    attr(pc.obj, "xCoord"), attr(pc.obj, "yCoord")
                              )
                        }else{
                              
                              pc.obj$Data[,,] <- mat2Dto3Darray(
                                    pc.obj[[1]][[1]]$orig,
                                    attr(pc.obj, "xCoord"), attr(pc.obj, "yCoord")
                              )
                        }
                  
                  attr(pc.obj$Data, "dimensions") = c("time", "lat", "lon")

            }
            
            
      pc.obj$xyCoords <- list( "x" = attr(pc.obj, "xCoord"), "y" = attr(pc.obj, "yCoord"))
      pc.obj$Dates[[1]]$start <- attr(pc.obj, "dates_start")
      pc.obj <- pc.obj[-(1:(nvars+1))]
      
      }else if(is.null(EOFs) | is.null(PCs)){
            stop("EOFs and/or PCs need to be defined")
      }else{
           pc.obj <- t(EOFs %*% t(PCs))     
      }
           return(pc.obj)
}


#end

