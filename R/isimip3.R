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
#' @keywords internal
#' @importFrom stats approxfun ecdf quantile
#' @impoFrom lubridate year
#' @author S. Herrera and M. Iturbide

isimip3 <- function(o, p, s, 
                     dates, 
                     lower_bound= c(NULL), 
                     lower_threshold= c(NULL), 
                     upper_bound= c(NULL), 
                     upper_threshold= c(NULL), 
                     randomization_seed= NULL, 
                     detrend= array(data = FALSE, dim = 1), 
                     rotation_matrices= c(NULL), 
                     n_quantiles=50, 
                     distribution= c("normal"), 
                     trend_preservation = array(data = "additive", dim=1), 
                     adjust_p_values = array(data = FALSE, dim=1), 
                     if_all_invalid_use = c(NULL), 
                     invalid_value_warnings=FALSE){
  meses.o <- months(as.Date(dates$obs_hist))
  meses.p <- months(as.Date(dates$sim_hist))
  meses.s <- months(as.Date(dates$sim_fut))
  meses.name <- unique(meses.s)
  years.o <- year(as.Date(dates$obs_hist))
  years.p <- year(as.Date(dates$sim_hist))
  years.s <- year(as.Date(dates$sim_fut))
  
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
    l <- adjust_bias_one_month(data, years, lower_bound = lower_bound, lower_threshold= lower_threshold, upper_bound= upper_bound, upper_threshold= upper_threshold,
                               randomization_seed = randomization_seed, detrend=detrend, rotation_matrices= rotation_matrices, n_quantiles=as.integer(n_quantiles), distribution= distribution,
                               trend_preservation =trend_preservation, adjust_p_values=adjust_p_values, if_all_invalid_use= if_all_invalid_use, invalid_value_warnings=invalid_value_warnings)
  })
  smap <- s
  for (m in c(1:length(meses.name))) {
    indMonth.s <- which(meses.s == meses.name[m])
    smap[indMonth.s] <- auxMonths[[m]][[1]]
  }
  return(smap)
}
#end
