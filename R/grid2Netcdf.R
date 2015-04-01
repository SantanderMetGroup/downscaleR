#' @title Write NetCDF
#' @description Export downscaleR object to NetCDF
#' 
#' @importFrom ncdf4 ncdim_def 
#' @importFrom ncdf4 ncvar_def 
#' @importFrom ncdf4 nc_create 
#' @importFrom ncdf4 ncatt_put 
#' @importFrom ncdf4 ncvar_put
#' @importFrom ncdf4 nc_close
#' 
#' @param data A grid data object coming from \code{\link{loadGridData}} or \code{\link{interpGridData}} 
#'  or the function \code{\link[ecomsUDG.Raccess]{loadECOMS}} of package \pkg{ecomsUDG.Raccess}. 
#' @param NetCDFOutFile Name of the file created by the function. ("out.nc4", default)
#' @param missval Value associated with the NA value used in R. (1e20, default)
#' @param globalAttributes List of global attributes included in the NetCDF file. (NULL, default)
#' @param varAttributes List of attributes to be included in the variable written in the NetCDF file. (NULL, default)
#' @param prec Precision to write the attribute. If not specified, the written precision is the same as the variable whose attribute this is. This can be overridden by specifying this argument with a value of "short", "float", "double", or "text".
#' @param compression If set to an integer between 1 (least compression) and 9 (most compression), this enables compression for the variable as it is written to the file. Turning compression on forces the created file to be in netcdf version 4 format, which will not be compatible with older software that only reads netcdf version 3 files. (4, default)
#' @param shuffle Turns on (if TRUE) or off (if FALSE, the default) the shuffle filter. According to netcdf docs, turning the shuffle filter on can improve compression for integer variables. Turning the shuffle filter on forces the created file to be in netcdf version 4 format, which will not be compatible with older software that only reads netcdf version 3 files.
#'  
#' @return A NetCDF-4 file with the variable and attributes defined in the inputs.
#' 
#' 
#' @references
#' 
#' \itemize{
#' \item David Pierce \email{dpierce@@ucsd.edu}, Interface to Unidata netCDF (version 4 or earlier) format data files, http://dwpierce.com/software
#' }
#' @author Sixto Herrera \email{herreras@@unican.es}, Wietse Franssen \email{wietse.franssen@@wur.nl}
#' @export
#' @examples \dontrun{
#' # The following panels show an illustrative use of ECOMS-UDG and downscaleR to obtain the bias corrected series of mean temperature for the period DJFMAM (one-month lead time; i.e. with the initializations from November) over Europe. WFDEI is used as reference. 
#' # Note that, in order to facilitate the use of the resulting bias corrected data in different impact applications, the resulting bias corrected data can be easily exported to netcdf format. 
#' 
#' obs <- loadECOMS(dataset = "WFDEI", var = "tas", season = c(12,1,2,3,4,5), lonLim = c(-15,35), latLim = c(32, 75), years = c(2001:2010))
#' prd <- loadECOMS(dataset = "System4_seasonal_15", var = "tas", time = "DD", season = c(12,1,2,3,4,5), members = 1:2, leadMonth = 1, lonLim = c(-15,35), latLim = c(32, 75), years = c(2001:2010))
#' prd <- interpGridData(prd, new.grid = getGrid(obs), method = "nearest");
#' 
#' # Bias correction and plotting
#' prd.bc <- biasCorrection(obs, prd, prd, method = "qqmap", multi.member = TRUE, window = 30) 
#' plotMeanField(obs) 
#' plotMeanField(prd, multi.member = FALSE) 
#' plotMeanField(prd.bc, multi.member = FALSE) 
#' 
#' # Exporting to netcdf4
#' fileName <- "tas_qqmap_System4_WFDEI_2001_2010.nc4"
#' grid2NetCDF(prd.bc, NetCDFOutFile = fileName, missval = 1e20, prec = "float")
#' # Including attributes:
#' globalAttributes <- list(institution = "Santander Meteorology Group, http://www.meteo.unican.es/")
#' varAttributes <- list(description = "Bias corrected daily mean temperature.")
#' grid2NetCDF(prd.bc, NetCDFOutFile = fileName, missval = 1e20, prec = "float", globalAttributes = globalAttributes, varAttributes = varAttributes)
#' }

grid2NetCDF <- function(data, NetCDFOutFile = "out.nc4", missval = 1e20, globalAttributes = NULL, varAttributes = NULL, prec = "float", compression=4, shuffle = TRUE) {
      data(vocabulary)
      ntime<-length(data$Dates$start)
      tmpStdName<-as.character(vocabulary$standard_name[pmatch(vocabulary$identifier,data$Variable$varName, nomatch=0) > 0])
      tmpUnits<-as.character(vocabulary$units[pmatch(vocabulary$identifier,data$Variable$varName, nomatch=0) > 0])
      time.index <- grep("^time$", attr(data$Data, "dimensions"))
      lon.index <- grep("^lon$", attr(data$Data, "dimensions"))
      lat.index <- grep("^lat$", attr(data$Data, "dimensions"))
      member.index <- grep("^member$", attr(data$Data, "dimensions"))
      datesList <- as.POSIXct(data$Dates$start, tz="GMT", format="%Y-%m-%d %H:%M:%S")
      times <- (as.double(datesList)-as.double(datesList[1]))/86400
      dimtime <- ncdim_def( "time", paste0("days since ", data$Dates$start[1]), times, unlim=FALSE, calendar="gregorian", create_dimvar = TRUE)
      dimlon  <- ncdim_def( "lon", units="degrees_east", data$xyCoords$x, longname="longitude", create_dimvar = TRUE)
      dimlat  <- ncdim_def( "lat", units="degrees_north", data$xyCoords$y, longname="latitude", create_dimvar = TRUE)
      if (length(member.index)>0){
            dimens  <- ncdim_def( "member", units="member", 0:(dim(data$Data)[member.index]-1), longname="realization", create_dimvar = TRUE)
            perOrdered <- c(lon.index,lat.index,member.index,time.index)
            dimOrdered <- list(dimlon,dimlat,dimens,dimtime)
      }else{
            perOrdered <- c(lon.index,lat.index,time.index)
            dimOrdered <- list(dimlon,dimlat,dimtime)
      }
      dataOrdered <- aperm(data$Data, perOrdered)
      var <- ncvar_def(data$Variable$varName, units=tmpUnits, dim=dimOrdered, missval, longname=tmpStdName, compression=compression, shuffle = shuffle)
      ncnew <- nc_create(NetCDFOutFile, var )
      ncatt_put(ncnew, "time", "standard_name","time")
      ncatt_put(ncnew, "time", "axis","T")
      ncatt_put(ncnew, "time", "_CoordinateAxisType","Time")
      ncatt_put(ncnew, "time", "_ChunkSize",1)
      ncatt_put(ncnew, "lon", "standard_name","longitude")
      ncatt_put(ncnew, "lon", "_CoordinateAxisType","Lon")
      ncatt_put(ncnew, "lat", "standard_name","latitude")
      ncatt_put(ncnew, "lat", "_CoordinateAxisType","Lat")
      if (length(member.index)>0){
            ncatt_put(ncnew, "member", "standard_name","realization")
            ncatt_put(ncnew, "member", "_CoordinateAxisType","Ensemble")
            ncatt_put(ncnew, "member", "ref","http://www.uncertml.org/samples/realisation")
      }
      ncatt_put(ncnew, data$Variable$varName, "missing_value",missval)
      for (v in 1:length(varAttributes)){
            ncatt_put(ncnew, data$Variable$varName, names(varAttributes)[v],as.character(varAttributes[v]))
      }
      for (v in 1:length(globalAttributes)){
            ncatt_put(ncnew, 0, names(globalAttributes)[v],as.character(globalAttributes[v]))
      }
      if (length(attr(data$Data, "correction"))>0){
            ncatt_put(ncnew, 0, "product","Bias-Correction")
            ncatt_put(ncnew, 0, "bc_method",attr(data$Data, "correction"))
      }
      if (length(attr(data, "dataset"))>0){
            ncatt_put(ncnew, 0, "dataset",attr(data, "dataset"))
      }
      if (length(attr(data, "source"))>0){
            ncatt_put(ncnew, 0, "source",attr(data, "source"))
      }
      ncatt_put(ncnew, 0, "description","NetCDF file created by downscaleR: https://github.com/SantanderMetGroup/downscaleR")
      ncatt_put(ncnew, 0, "Conventions","CF-1.4")
      ncvar_put(ncnew, var, dataOrdered )
      nc_close(ncnew)
      message(paste0("NetCDF file written: ", NetCDFOutFile))
}
