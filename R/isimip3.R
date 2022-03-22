#' @title ISIMIP v3 method for bias correction
#' @description Implementation of ISIMIP v3 method for bias correction:
#'   1. Replaces invalid values in time series.
#'   2. Detrends time series if desired.
#'   3. Replaces values beyond thresholds by random numbers.
#'   4. Adjusts inter-variable copula.
#'   5. Adjusts marginal distributions for every variable.
#'   6. Replaces values beyond thresholds by the respective bound.
#'   7. Restores trends.
#'   Removes invalid values from masked arrays and store resulting numpy arrays
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param s A vector containing the simulated climate for the variable used in \code{p}, but considering the test period.
#' @param dates List containing the dates corresponding to Observations (obs_hist), training period (sim_hist) and test period (sim_fut). 
#'        Every array represents the years of the time steps of the time series data. Used for detrending.
#' @param long_term_mean : dict of lists with keys 'obs_hist', 'sim_hist', 'sim_fut'
#'        Every value in every list represents the average of all valid values in
#'        the complete time series for one climate variable.
#' @param lower_bound : list of floats, optional. Lower bounds of values in x_obs_hist, x_sim_hist, and x_sim_fut.
#' @param lower_threshold : list of floats, optional. Lower thresholds of values in x_obs_hist, x_sim_hist, and x_sim_fut. 
#'        All values below this threshold are replaced by random numbers between lower_bound and lower_threshold before bias adjustment.
#' @param upper_bound : list of floats, optional. Upper bounds of values in x_obs_hist, x_sim_hist, and x_sim_fut.
#' @param upper_threshold : list of floats, optional. Upper thresholds of values in x_obs_hist, x_sim_hist, and x_sim_fut. 
#'        All values above this threshold are replaced by random numbers between upper_threshold and upper_bound before bias adjustment.
#' @param randomization_seed : int, optional. Used to seed the random number generator before replacing invalid values and values beyond 
#'        the specified thresholds.
#' @param detrend : list of booleans, optional. Detrend time series before bias adjustment and put trend back in afterwards.
#' @param rotation_matrices : list of (n,n) ndarrays, optional. List of orthogonal matrices defining a sequence of rotations in variable 
#'        space, where n is the number of variables.
#' @param n_quantiles : int, optional. Number of quantile-quantile pairs used for non-parametric quantile mapping.
#' @param distribution : list of strs, optional. Kind of distribution used for parametric quantile mapping:
#'        {"normal", "weibull", "gamma", "beta", "rice"}.
#' @param trend_preservation : list of strs, optional. Kind of trend preservation used for non-parametric quantile mapping: 
#'        {"additive", "multiplicative", "mixed", "bounded"}.
#' @param adjust_p_values : list of booleans, optional. Adjust p-values for a perfect match in the reference period.
#' @param if_all_invalid_use : list of floats, optional. Used to replace invalid values if there are no valid values. An error is raised 
#'        if there are no valid values and this parameter is None.
#' @param invalid_value_warnings : boolean, optional. Raise user warnings when invalid values are replaced bafore bias adjustment.
#' @param halfwin_upper_bound_climatology : Determines the lengths of running windows used in the calculations of climatologies of upper bounds that are used to scale values 
#' of obs_hist, sim_hist, and sim_fut to the interval [0,1] before bias adjustment. The window length is set to: 
#' halfwin_upper_bound_climatology * 2 + 1 time steps. If halfwin_upper_bound_climatology == 0 then no rescaling is done.
#' @keywords internal
#' @importFrom stats approxfun ecdf quantile
#' @importFrom reticulate dict source_python
#' @author S. Herrera and M. Iturbide

isimip3 <- function(o, p, s,
                    dates, 
                    long_term_mean = NULL,
                    lower_bound = c(NULL),
                    lower_threshold = c(NULL),
                    upper_bound = c(NULL),
                    upper_threshold = c(NULL),
                    randomization_seed = c(NULL),
                    detrend= array(data = FALSE, dim = 1),
                    rotation_matrices = c(NULL),
                    n_quantiles = 50,
                    distribution = c("normal"),
                    trend_preservation = array(data = "additive", dim = 1),
                    adjust_p_values = array(data = FALSE, dim = 1),
                    if_all_invalid_use = c(NULL),
                    invalid_value_warnings = FALSE,
                    halfwin_upper_bound_climatology = 0) {
      
      smap <- s
      if (any(!is.na(o)) & any(!is.na(p)) & any(!is.na(s))) {
             ## source python routines
            lf <- list.files(file.path(find.package("downscaleR")), pattern = "\\.py$", recursive = TRUE, full.names = TRUE)
            sapply(lf, source_python, .GlobalEnv)
            if (halfwin_upper_bound_climatology != 0){
                  aux <- scale_by_upper_bound_climatology(o,dates$obs_hist, halfwin = halfwin_upper_bound_climatology)
                  o <- aux$data
                  aux <- scale_by_upper_bound_climatology(p,dates$sim_hist, halfwin = halfwin_upper_bound_climatology)
                  p <- aux$data
                  s.rescale <- scale_by_upper_bound_climatology(s,dates$sim_fut, halfwin = halfwin_upper_bound_climatology)
                  s <- s.rescale$data
            }
            if (is.null(long_term_mean)) {
                  lt.o <- average_valid_values(array(o), if_all_invalid_use = if_all_invalid_use,
                                               lower_bound=lower_bound, lower_threshold = lower_threshold,
                                               upper_bound=upper_bound, upper_threshold = upper_threshold)
                  lt.p <- average_valid_values(array(p), if_all_invalid_use = if_all_invalid_use,
                                               lower_bound=lower_bound, lower_threshold = lower_threshold,
                                               upper_bound=upper_bound, upper_threshold = upper_threshold)
                  lt.s <- average_valid_values(array(s), if_all_invalid_use = if_all_invalid_use,
                                               lower_bound=lower_bound, lower_threshold = lower_threshold,
                                               upper_bound=upper_bound, upper_threshold = upper_threshold)
                  long_term_mean <- list(obs_hist = lt.o, sim_hist = lt.p, sim_fut = lt.s)
            }
            meses.o <- as.numeric(substr(dates$obs_hist, 6, 7))
            meses.p <- as.numeric(substr(dates$sim_hist, 6, 7))
            meses.s <- as.numeric(substr(dates$sim_fut, 6, 7))
            meses.name <- unique(meses.s)
            years.o <- as.numeric(substr(dates$obs_hist, 1, 4))
            years.p <- as.numeric(substr(dates$sim_hist, 1, 4))
            years.s <- as.numeric(substr(dates$sim_fut, 1, 4))
            
            ## Loop in months:
            auxMonths <- lapply(1:length(meses.name), function(m) {
                  indMonth.o <- which(meses.o == meses.name[m])
                  indMonth.p <- which(meses.p == meses.name[m])
                  indMonth.s <- which(meses.s == meses.name[m])
                  data <- dict(obs_hist = matrix(o[indMonth.o], ncol = length(indMonth.o), nrow = 1),
                               sim_hist = matrix(p[indMonth.p], ncol = length(indMonth.p), nrow = 1),
                               sim_fut = matrix(s[indMonth.s], ncol = length(indMonth.s), nrow = 1))
                  years <- dict(obs_hist = matrix(years.o[indMonth.o], ncol = length(indMonth.o), nrow = 1),
                                sim_hist = matrix(years.p[indMonth.p], ncol = length(indMonth.p), nrow = 1),
                                sim_fut = matrix(years.s[indMonth.s], ncol = length(indMonth.s), nrow = 1))
                  l <- adjust_bias_one_month(data, years, long_term_mean, lower_bound = lower_bound, lower_threshold= lower_threshold, upper_bound= upper_bound, upper_threshold= upper_threshold,
                                             randomization_seed = randomization_seed, detrend=detrend, rotation_matrices= rotation_matrices, n_quantiles=as.integer(n_quantiles), distribution= distribution,
                                             trend_preservation =trend_preservation, adjust_p_values=adjust_p_values, invalid_value_warnings=invalid_value_warnings)
            })
            for (m in c(1:length(meses.name))) {
                  indMonth.s <- which(meses.s == meses.name[m])
                  smap[indMonth.s] <- auxMonths[[m]][[1]]
            }
            if (halfwin_upper_bound_climatology != 0){
                  smap <- re_scale_by_upper_bound_climatology(smap,dates$sim_fut,s.rescale)
            }
      } else{
             smap[] <- NA
##             print("No valid values have been found, so the raw un-corrected data is returned.") 
      }
      return(smap)
}
#end

#' @title Scaling by upper bound climatology
#' @description Scaling by upper bound climatology
#' @param data A vector containing the climate data
#' @param dates A vector containing the dates to which corresponds the data. 
#' @param halfwin A number defining the length of running windows used in the calculations of 
#'                climatologies of upper bounds that are used to scale values. The
#'                window length is set to halfwin * 2 + 1 time
#'                steps. If halfwin == 0 then no rescaling is done.
#' 
#' @keywords internal
#' @author S. Herrera (sixto.herrera@@unican.es)

scale_by_upper_bound_climatology <- function(data, dates, halfwin = 0) {
      if (halfwin != 0){
        doys.o <- unique(substr(dates, 6, 10))
        mydmrm.o <- lapply(c(1:length(doys.o)), function(d){
          ind.d <- which(substr(dates, 6, 10) == doys.o[d])
          if (length(ind.d) > 0){
            d.p <- lapply(c(1:length(ind.d)), function(dd){
              ind.dd <- c(max(c(1,ind.d[dd]-halfwin)):min(c(length(data),ind.d[dd]+halfwin)))
              l <- max(data[ind.dd], na.rm = TRUE)
            })
            l <- max(as.numeric(unlist(d.p)), na.rm = TRUE)
          }else{
            l <- NA
          }
        })
        mydmrm.o <- as.numeric(unlist(mydmrm.o))
        ubc.o <- mydmrm.o
        mydmrm.o <- c(mydmrm.o[(length(doys.o)-halfwin+1):length(doys.o)],mydmrm.o,mydmrm.o[1:halfwin])
        for (d in c(1:length(doys.o))){
          center <- d+halfwin
          ubc.o[d] <- mean(mydmrm.o[c((center-halfwin):(center+halfwin))], na.rm = TRUE)
          ind.d <- which(substr(dates, 6, 10) == doys.o[d])
          if (length(ind.d) > 0){
            data[ind.d] <- data[ind.d]/ubc.o[d]
          }
        }
      }
      out <- list("data" = data, "doys" = doys.o, "normd" = ubc.o)
      return(out)
    }
#end

#' @title Re-Scaling by upper bound climatology
#' @description Re-Scaling by upper bound climatology
#' @param data A vector containing the climate data
#' @param dates A vector containing the dates to which corresponds the data. 
#' @param rescale A list with the elements needed (day-of-year and scale corresponding to each day-of-year) to re-scale the bias adjusted data.
#' 
#' @keywords internal
#' @author S. Herrera (sixto.herrera@@unican.es)

re_scale_by_upper_bound_climatology <- function(data, dates, rescale) {
      for (d in c(1:length(rescale$doys))){
        ind.d <- which(substr(dates, 6, 10) == rescale$doys[d])
        if (length(ind.d) > 0){
          data[ind.d] <- data[ind.d]*rescale$normd[d]
        }
      }
      return(data)
    }
#end
