downscale.chunk <- function(x, y, newdata,
                            method, ...,
                            global.vars = NULL, combined.only = TRUE, spatial.predictors = NULL, local.predictors = NULL, extended.predictors = NULL,
                            condition = NULL, threshold = NULL,
                            path = getwd()) {
  
  x <- getTemporalIntersection(x,y,which.return = "obs")
  y <- getTemporalIntersection(x,y,which.return = "prd")
  
  latitudes <- getCoordinates(y)$y
  chunks <- length(latitudes)
  lapply(1:chunks,FUN = function(z){
    print(paste("Training chunk:",z,"out of",chunks))
    y_chunk <- subsetDimension(y,dimension = "lat", indices = z)
    xyT <- prepareData(x = x, y = y_chunk, global.vars = global.vars, combined.only = combined.only, spatial.predictors = spatial.predictors, local.predictors = local.predictors, extended.predictors = extended.predictors)
    model <- downscale.train(xyT, method, condition, threshold, ...)
    
    p <- lapply(newdata, function(zz) {
      xyt <- prepareNewData(zz,xyT)
      downscale.predict(newdata = xyt,model)
    })
    
    if (z < 10) {zn <- paste0(0,z)}
    else {zn <- z}
    lapply(1:(length(p)+1), function(zzz) {
    # lapply(2:(length(p)+1), function(zzz) {
    # lapply(1:1, function(zzz) {
      if (zzz == 1) {
        grid <- model$pred
      }
      else{
        grid <- p[[zzz-1]] 
      }
      save(grid, file = paste0(path,"/dataset",zzz,"_chunk",zn,".rda"))
    })
    p <- NULL
  })
  NULL
}
