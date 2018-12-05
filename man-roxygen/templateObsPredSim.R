#' @param y A grid or station data containing the observed climate data for the training period
#' @param x A grid containing the simulated climate by the model for the training period. This can be either
#'  the same variable as \code{obs}, in the case of model calibration (bias correction and related techniques) or MOS (model
#'   output statistics) approaches, or a set of predictors in the case of \emph{perfect prog} downscaling approaches (possibly
#'   after principal component analysis via the \code{\link{prinComp}} output).
#' @param newdata A grid containing the simulated climate for the variables used in \code{x}, but considering the test period.
#' If \code{x} is a multimember grid (the member dimension is > 1), the function assumes that the order of 
#' the members is the same in \code{newdata}.

