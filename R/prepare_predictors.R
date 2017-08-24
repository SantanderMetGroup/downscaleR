#   prepare_predictors.R Configuration of predictors for downscaling
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

#' @title Configuration of predictors for downscaling
#' @description Configuration of predictors for flexible downscaling experiment definition
#' @param x A grid (usually a multigrid) of predictor fields
#' @param subset.vars An optional character vector with the short names of the variables of the input \code{x} 
#'  multigrid to be retained as global predictors (use the \code{\link{getVarNames}} helper if not sure about variable names).
#'  This argument just produces a call to \code{\link[transformeR]{subsetGrid}}, but it is included here for better
#'  flexibility in downscaling experiments (predictor screening...). For instance, it allows to use some 
#'  specific variables contained in \code{x} as local predictors and the remaining ones, specified in \code{subset.vars},
#'  as either raw global predictors or to construct the combined PC.
#' @param y A grid (usually a stations grid, but not necessarily) of observations (predictands)
#' @param PCA Default to \code{NULL}, and not used. Otherwise, a named list of arguments in the form \code{argument = value},
#'  with the arguments to be passed to \code{\link[transformeR]{prinComp}} to perform Principal Component Analysis
#'  of the predictors grid (\code{x}). See Details on principal component analysis of predictors.
#' @param local.predictors Default to \code{NULL}, and not used. Otherwise, a named list of arguments in the form \code{argument = value},
#'  with the following arguments:
#'  \itemize{
#'    \item \code{neigh.vars}: names of the variables in \code{x} to be used as local predictors
#'    \item \code{neigh.fun}: Optional. Aggregation function for the selected local neighbours.
#'    The aggregation function is specified as a list, indicating the name of the aggregation function in
#'     first place (as character), and other optional arguments to be passed to the aggregation function.
#'     For instance, to compute the average skipping missing values: \code{neigh.fun = list(FUN= "mean", na.rm = TRUE)}.
#'     Default to NULL, meaning that no aggregation is performed.
#'    \item \code{n.neighs}: Integer. Number of nearest neighbours to use. If a single value is introduced, and there is more
#'    than one variable in \code{neigh.vars}, the same value is used for all variables. Otherwise, this should be a vector of the same
#'    length as \code{neigh.vars} to indicate a different number of nearest neighbours for different variables.
#'  }
#'  
#' @return A named list with components \code{y} (the predictand), \code{x.global} (global predictors, 2D matrix), \code{x.local} (local predictors, a list) 
#' and \code{pca} (\code{\link[transformeR]{prinComp}} output), and other attributes. See Examples.
#'  
#' @examples 
#' # See the dedicated vignette by typing:
#' # utils::vignette(topic = "configuring_predictors", package = "downscaleR")
#'  
#' @details   
#'  \strong{Temporal consistency}
#'  Note that \code{x} (predictors) and \code{y} predictands are checked for temporal consistency
#'   prior to downscaling. In case of partial temporal overlapping, both are internnaly intersected for exact temporal matching.
#'  
#'  \strong{Principal Component Analysis}
#'  Always that PCA is used, a combined PC will be returned (unless one single predictor is used, case in which no combination is possible).
#'  Note that the variables of the predictor grid used to construct the combined PC can be flexibly controlled through the optional argument
#'  \code{subset.vars}.
#'  
#' @importFrom transformeR getTemporalIntersection getRefDates getCoordinates getVarNames
#' @importFrom magrittr %<>% %>% 
#'  
#' @family downscaling.helpers
#'  
#' @export
#'  
#' @author J. Bedia, D. San-Mart\'in and J.M. Guti\'errez 

prepare_predictors <- function(x, y, subset.vars = NULL, PCA = NULL, local.predictors = NULL) {
    y <- getTemporalIntersection(obs = y, prd = x, which.return = "obs")
    x <- getTemporalIntersection(obs = y, prd = x, which.return = "prd")
    dates <- getRefDates(x)
    xyCoords <- getCoordinates(x)
    # Local predictors
    nn <- NULL
    if (!is.null(local.predictors)) {
            pred.nn.ind <- predictor.nn.indices(neigh.vars = local.predictors$neigh.vars,
                                                n.neighs = local.predictors$n.neighs,
                                                x = x, y = y) 
            nn <- predictor.nn.values(nn.indices.list = pred.nn.ind,
                                      grid = x,
                                      neigh.fun = local.predictors$neigh.fun) 
    } 
    # Global predictor subsetting
    if (!is.null(subset.vars)) {
        x  %<>%  subsetGrid(var = subset.vars, drop = FALSE)
    }
    # PCA - only considers the combined PC as global predictor
    pc.comb <- NULL
    full.pca <- NULL
    if (!is.null(PCA)) {
        if ("which.var" %in% names(PCA)) {
            PCA[["which.var"]] <- subset.vars
            message("NOTE: The argument 'which.var' passed to the 'PCA' list was overriden by argument 'subset.vars'")
        }
        varnames <- getVarNames(x)
        PCA[["grid"]] <- x
        PCA[["combined.PC"]] <- TRUE
        full.pca <- do.call("prinComp", PCA)
        if (length(varnames) == 1) {
            pc.comb <- full.pca %>% extract2(varnames)
        } else {
            pc.comb <- full.pca %>% extract2("COMBINED")    
        }
        x <- pc.comb[[1]]$PCs
    # RAW predictors    
    } else {
        # Construct a predictor matrix of raw input grids
        nvars <- getShape(x, "var")
        varnames <- getVarNames(x)
        aux.list <- lapply(1:nvars, function(i) {
            subsetGrid(x, var = varnames[i]) %>% redim(drop = TRUE) %>% extract2("Data") %>% array3Dto2Dmat()
        })
        x <- do.call("cbind", aux.list)
        aux.list <- NULL
    }
    predictor.list <- list("y" = y, "x.global" = x, "x.local" = nn, "pca" = full.pca)
    if (!is.null(PCA)) PCA <- PCA[-grep("grid", names(PCA))]
    attr(predictor.list, "PCA.pars") <- PCA
    attr(predictor.list, "localPred.pars") <- local.predictors
    attr(predictor.list, "dates") <- dates
    attr(predictor.list, "xyCoords") <- xyCoords
    return(predictor.list)
}


#' @title Calculate spatial index position of nearest neighbors
#' @description Calculate spatial index position of nearest neighbors
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
    }
    if (length(n.neighs) != length(neigh.vars)) {
        stop("Incorrect number of neighbours selected: this should be either 1 or a vector of the same length as 'neigh.vars'", call. = FALSE)
    }
    # The index list has the same length as the number of local predictor variables, containing a matrix of index positions
    local.pred.list <- lapply(1:length(neigh.vars), function(j) {
        ind.mat <- vapply(1:nrow(coords.y), FUN.VALUE = numeric(n.neighs[j]), FUN = function(i) {
            dists <- sqrt((coords.y[i,1] - coords.x[,1])^2 + (coords.y[i,2] - coords.x[,2])^2)
            which(dists %in% sort(dists)[1:n.neighs[j]])
        })
    })
    names(local.pred.list) <- neigh.vars
    return(local.pred.list)
}


#' @title Construct local predictor matrices
#' @description Constructs the local predictor matrices given their spatial index position 
#' (as returned by \code{\link{predictor.nn.indices}}) and additional parameters.
#' @param nn.indices.list A list of index positions as returned by \code{\link{predictor.nn.indices}}
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
    } 
    if (!is.null(neigh.fun) && (length(neigh.fun) != length(nn.indices.list))) {
        stop("Incorrect number of neighbours selected: this should be either 1 or a vector of the same length as 'nn.indices.list'", call. = FALSE)
    }
    out.list <- lapply(1:length(nn.indices.list), function(j) {
        neigh.vars <- names(nn.indices.list)
        aux <- subsetGrid(grid, var = neigh.vars[j], drop = TRUE) %>% extract2("Data") %>% array3Dto2Dmat()
        # Local predictor matrix list
        l <- lapply(1:ncol(nn.indices.list[[j]]), function(k) aux[ ,nn.indices.list[[j]][ ,k]])
        # Neighbour aggregation function 
        if (!is.null(neigh.fun[[j]])) {
            arg.list <- c(neigh.fun[[j]], "MARGIN" = 1)
            lapply(1:length(l), function(k) do.call(what = "apply", c("X" = l[k], arg.list))) 
        } else {
            l
        }
    })
    lapply(1:length(out.list[[1]]), function(x) {
        expr <- paste0("cbind(", paste0("out.list[[", 1:length(out.list), "]][[x]]", collapse = ","), ")")
        parse(text = expr) %>% eval()
    })
}


