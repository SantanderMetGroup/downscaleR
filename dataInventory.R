# dataset can be either a ncml-like file or a directory in the case of observations
dataInventory <- function(dataset, return.stats = FALSE) {
      rs <- return.stats
      message(paste("[", Sys.time(), "] Doing inventory ...", sep = ""))
      if (isTRUE(file.info(dataset)$isdir)) {
            out <- dataInventory.ASCII(dataset, rs)
      } else {
            out <- dataInventory.NetCDF(dataset)
      }
      message(paste("[", Sys.time(), "] Done.", sep = ""))
      return(out)
}
# End


# # Examples
# # 1. NCEP reanalysis
# dataset = '/home/juaco/workspace/gitRepo/downscaling/datasets/reanalysis/Iberia_NCEP/Iberia_NCEP.ncml'
# station.inv <- dataInventory(dataset)
# str(station.inv)
# names(station.inv)
# station.inv$Q$Description

# # 2. Seasonal Forecast
# dataset = "http://www.meteo.unican.es/tds5/dodsC/system4/System4_Seasonal_51Members.ncml"
# source("../useRportal/loginTHREDDS.R")
# loginTHREDDS("juaco", "********")
# sys4.inventory <- dataInventory(dataset) # Note that a summary table is printed on screen
# # Names of the variables
# names(sys4.inventory) 
# # Structure of the inventory
# str(sys4.inventory) 
# # Details of variable sea level pressure
# # Initializations (Run times)
# sys4.inventory$Mean_sea_level_pressure_surface$Dimensions$run
# # Forecast Time range
# range(sys4.inventory$Mean_sea_level_pressure_surface$Dimensions$time1$Values)


