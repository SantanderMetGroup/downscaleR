# data("NCEP_Iberia_hus850", "NCEP_Iberia_psl", "NCEP_Iberia_ta850", "NCEP_Iberia_tas")
# x <- makeMultiGrid(NCEP_Iberia_hus850, NCEP_Iberia_psl, NCEP_Iberia_ta850)
# data("VALUE_Iberia_tp")
# y <- VALUE_Iberia_tp
# newdata = NULL
# 
# getVarNames(x)
# 
# 
# str(x)
# 
# "PCA" = list(n.eofs = c(5,5,3,8),
#              v.exp = .975,
#              combined.PC = TRUE,
#              which.combine = c("hus850", "psl")),
# "local.predictors" = list(neigh.vars = "hus850",
#                           n.neighs = 5,
#                           neigh.fun = NULL)


#   prepare.predictors.R Configuration of predictors for downscaling
#
#   Copyright (C) 2017 Santander Meteorology Group (http://www.meteo.unican.es)
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
# 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
# 
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

#'  @title Configuration of predictors for downscaling
#'  @description Configuration of predictors for flexible downscaling experiment definition
#'  @param x A grid (usually a multigrid) of predictor fields
#'  @param y A grid (usually a stations grid, but not necessarily) of observations (predictands)
#'  @param PCA Default to \cde{NULL}, and not used. Otherwise, a named list of arguments in the form \code{argument = value},
#'  with the arguments to be passed to \code{\link[transformeR]{prinComp}} to perform Principal Component Analysis
#'  of the predictors field. See Details on principal component analysis of predictors
#'  @param local.predictors Default to \cde{NULL}, and not used. Otherwise, a named list of arguments in the form \code{argument = value},
#'  with the arguments
#'  
#'  @details   
#'  \strong{Temporal consistency}
#'  Note that \code{x} (predictors) and \code{y} predictands are checked for temporal consistency
#'   prior to downscaling. In case of partial temporal overlapping, both are internnaly intersected for exact temporal matching.
#'  
#'  \strong{Principal Component Analysis}
#'  Always that PCA is used, a combined PC will be returned (unless one single predictor is used, case in which no combination is possible).
#'  Note that the variables of the predictor grid used to constrcut the combined PC can be flexibly controlled trough the optional argument
#'  \code{which.combine} in \code{\link[transformeR]{prinComp}}.
#'  
#'  @author J. Bedia, D. San-Mart\'in and J.M. Guti\'errez 
#'  @importFrom transformeR getTemporalIntersection
#'  @export

prepare.predictors <- function(x, y, PCA = NULL, local.predictors = NULL) {
    y <- getTemporalIntersection(obs = y, prd = x, which.return = "obs")
    x <- getTemporalIntersection(obs = y, prd = x, which.return = "prd")
    nn <- if (!is.null(local.predictors)) {
        pred.nn.ind <- predictor.nn.indices(neigh.vars = local.predictors$neigh.vars,
                                            n.neighs = local.predictors$n.neighs,
                                            x = x, y = y) 
        predictor.nn.values(nn.indices.list = pred.nn.ind,
                            grid = x,
                            neigh.fun = local.predictors$neigh.fun) 
    } else {
        NULL
    }
    if (!is.null(PCA)) {
        do.call("prinComp", PCA) %>% extract2("COMBINED")
    }
}


# load("ignore/juaco/data/obsPredSim_NCEP.Rdata", verbose = TRUE)
# # obs <- obs.tmean
# 
# data("VALUE_Iberia_tas")
# obs <- VALUE_Iberia_tas
# 
# getShape(pred)
# getShape(sim)

# Check the consistency of x and newdata
# checkDim(pred, sim, dimensions =  c("var","lon", "lat"))

# Obtain the coordinates of x (=newdata) and y
# coords.y <- get2DmatCoordinates(obs)
# coords.x <- get2DmatCoordinates(pred)

# Compute PCA
# xpc <- prinComp(pred, v.exp = .9)
# 
# # Target variable(s) of neighbors (either number or name)
# neigh.vars = c("psl", "tas")
# n.neigh = c(9,3) # number of neighbors
# neigh.fun = list(list(FUN = "mean", na.rm = TRUE), NULL)
# 
# a <- predictor.nn.indices(neigh.vars = c("psl", "tas"), n.neighs = c(10,3), x.grid = pred, y.grid = obs)
# q <- predictor.nn.values(a, neigh.fun = list(NULL, list(FUN = mean, na.rm = TRUE)), grid = pred)
# 
# str(q)

#' @param neigh.vars Character vector of the names of the variable(s) to be used as local predictors
#' @param n.neighs Integer vector. Number of nearest neighbours to be considered for each variables in \code{neigh.vars}.
#' Its length must be either 1 --in case there is one single local predictor variable or the same number of neighbours is 
#' desired for every local predictor variable-- or equal the length of \code{neigh.vars}, to consider a varying number of 
#' neighbours for each local predictor variable.
#' @param x.grid The (multi)grid containing the predictors
#' @param y.grid The grid containing the predictand
#' @keywords internal
#' @importFrom transformeR getVarNames 
#' @return A list of matrices (rows = local neighbors, cols = predictand points). The list is of the same
#' length as the number of local predictor variables chosen
#' @family downscaling.helpers
#' @author J Bedia

predictor.nn.indices <- function(neigh.vars = NULL, n.neighs = NULL, x, y) {
    if (is.null(neigh.vars)) stop("Undefined local predictor variables. A value for 'neigh.vars' argument is required", call. = FALSE)
    if (is.null(n.neighs)) stop("Undefined number of local neighbours. A value for 'n.neighs' argument is required", call. = FALSE)
    varnames <- getVarNames(x)
    if (!any(neigh.vars %in% varnames)) stop("The requested neigh.var was not found in the predictor grid", call. = FALSE)
    # Selection of nearest neighbour indices
    # Coordinates matrices
    coords.y <- get2DmatCoordinates(y)
    coords.x <- get2DmatCoordinates(x)
    if (any(n.neighs > nrow(coords.x))) stop("Too many neighbours selected (more than predictor grid cells)", call. = FALSE) 
    if (length(n.neighs) == 1) {
        n.neighs <- rep(n.neighs, length(neigh.vars))
    } else if (length(n.neighs) != length(neigh.vars)) {
        stop("Incorrect number of neighbours selected: this should be either 1 or a vector of the same length as 'neigh.vars'", call. = FALSE)
    }
    # The index list has the same length as the number of local predictor variables, containing a matrix of index positions
    local.pred.list <- lapply(1:length(neigh.vars), function(j) {
        ind.mat <- vapply(1:nrow(coords.y), numeric(n.neighs[j]), FUN = function(i) {
            dists <- sqrt((coords.y[i,1] - coords.x[,1])^2 + (coords.y[i,2] - coords.x[,2])^2)
            which(dists %in% sort(dists)[1:n.neighs[j]])
        })
    })
    names(local.pred.list) <- neigh.vars
    return(local.pred.list)
}



#' @param nn.indices.list
#' @param neigh.fun A list of funtions for neighbour aggregation for each local predictor variable.
#'  Default to \code{NULL}, and no aggregation is performed
#' @param grid The grid containing the local predictor information. It can be either the predictors (training)
#'  grid or the simulation (prediction) grid, depending on the situation.
#' @keywords internal
#' @family downscaling.helpers
#' @importFrom transformeR subsetGrid array3Dto2Dmat
#' @importFrom magrittr %>% extract2
#' @author J Bedia

predictor.nn.values <- function(nn.indices.list, grid, neigh.fun = NULL) {
    # Aggregation function of neighbors
    if (length(nn.indices.list) == 1) {
        neigh.fun <- list(neigh.fun)
    } else if (length(neigh.fun) != length(nn.indices.list)) {
        stop("Incorrect number of neighbours selected: this should be either 1 or a vector of the same length as 'nn.indices.list'", call. = FALSE)
    }
    out.list <- lapply(1:length(nn.indices.list), function(j) {
        neigh.vars <- names(nn.indices.list)
        aux <- subsetGrid(grid, var = neigh.vars[j], drop = TRUE) %>% extract2("Data") %>% array3Dto2Dmat()
        # Local predictor matrix list
        l <- lapply(1:ncol(nn.indices.list[[j]]), function(k) aux[,nn.indices.list[[j]][,k]])
        # Neighbour aggregation function 
        if (!is.null(neigh.fun[[j]])) {
            arg.list <- c(neigh.fun[[j]], "MARGIN" = 1)
            lapply(1:length(l), function(k) do.call(what = "apply", c("X" = l[k], arg.list))) 
        } else {
            l
        }
    })
    lapply(1:length(out.list[[1]]), function(x) cbind(out.list[[1]][[x]], out.list[[2]][[x]]))
}


