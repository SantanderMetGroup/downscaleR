#   prepareData.R Configuration of predictors for downscaling
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

#' @title Configuration of data for downscaling
#' @description Configuration of data for flexible downscaling experiment definition
#' @param x A grid (usually a multigrid) of predictor fields
#' @param global.vars An optional character vector with the short names of the variables of the input \code{x} 
#'  multigrid to be retained as global predictors (use the \code{\link{getVarNames}} helper if not sure about variable names).
#'  This argument just produces a call to \code{\link[transformeR]{subsetGrid}}, but it is included here for better
#'  flexibility in downscaling experiments (predictor screening...). For instance, it allows to use some 
#'  specific variables contained in \code{x} as local predictors and the remaining ones, specified in \code{subset.vars},
#'  as either raw global predictors or to construct the combined PC.
#' @param y A grid (usually a stations grid, but not necessarily) of observations (predictands)
#' @param spatial.predictors Default to \code{NULL}, and not used. Otherwise, a named list of arguments in the form \code{argument = value},
#'  with the arguments to be passed to \code{\link[transformeR]{prinComp}} to perform Principal Component Analysis
#'  of the predictors grid (\code{x}). See Details on principal component analysis of predictors.
#' @param combined.only Optional, and only used if spatial.predictors parameters are passed. Should the combined PC be used as the only
#' global predictor? Default to TRUE. Otherwise, the combined PC constructed with \code{which.combine} argument in 
#' \code{\link{prinComp}} is append to the PCs of the remaining variables within the grid.
#' @param local.predictors Default to \code{NULL}, and not used. Otherwise, a named list of arguments in the form \code{argument = value},
#'  with the following arguments:
#'  \itemize{
#'    \item \code{vars}: names of the variables in \code{x} to be used as local predictors
#'    \item \code{fun}: Optional. Aggregation function for the selected local neighbours.
#'    The aggregation function is specified as a list, indicating the name of the aggregation function in
#'     first place (as character), and other optional arguments to be passed to the aggregation function.
#'     For instance, to compute the average skipping missing values: \code{fun = list(FUN= "mean", na.rm = TRUE)}.
#'     Default to NULL, meaning that no aggregation is performed.
#'    \item \code{n}: Integer. Number of nearest neighbours to use. If a single value is introduced, and there is more
#'    than one variable in \code{vars}, the same value is used for all variables. Otherwise, this should be a vector of the same
#'    length as \code{vars} to indicate a different number of nearest neighbours for different variables.
#'  }
#' @param extended.predictors This is a parameter related to the extreme learning machine and reservoir computing framework where input data is randomly projected into a new space of size \code{n}. Default to \code{NULL}, and not used. Otherwise, a named list of arguments in the form \code{argument = value},
#'  with the following arguments:
#'  \itemize{
#'    \item \code{n}: A numeric value. Indicates the size of the random nonlinear dimension where the input data is projected.
#'    \item \code{module}: A numeric value (Optional). Indicates the size of the mask's module. Belongs to a specific type of ELM called RF-ELM.
#'  }
#' @return A named list with components \code{y} (the predictand), \code{x.global} (global predictors, 2D matrix), \code{x.local} (local predictors, a list) 
#' and \code{pca} (\code{\link[transformeR]{prinComp}} output), and other attributes. See Examples.
#'  
#' @details   
#'  \strong{Temporal consistency}
#'  Note that \code{x} (predictors) and \code{y} predictands are checked for temporal consistency
#'   prior to downscaling. In case of partial temporal overlapping, both are internnaly intersected for exact temporal matching.
#'  
#'  \strong{Principal Component Analysis}
#'  Always that spatial.predictors is used, a combined PC will be returned (unless one single predictor is used, case in which no combination is possible).
#'  Note that the variables of the predictor grid used to construct the combined PC can be flexibly controlled through the optional argument
#'  \code{subset.vars}.
#'  
#' @importFrom transformeR getTemporalIntersection getRefDates getCoordinates getVarNames
#' @importFrom magrittr %<>% %>% 
#' @seealso \href{https://github.com/SantanderMetGroup/downscaleR/wiki/preparing-predictor-data}{downscaleR Wiki} for preparing predictors for downscaling and seasonal forecasting.
#' @family downscaling.helpers
#' @export
#'  
#' @author J. Bedia, D. San-Martín and J.M. Gutiérrez 
#' @examples
#' # Loading data
#' data("VALUE_Iberia_tas")
#' y <- VALUE_Iberia_tas 
#' data("NCEP_Iberia_hus850", "NCEP_Iberia_psl", "NCEP_Iberia_ta850")
#' x <- makeMultiGrid(NCEP_Iberia_hus850, NCEP_Iberia_psl, NCEP_Iberia_ta850)
#' # Raw data
#' data <- prepareData(x = x, y = y)
#' # Using PCs as predictors. Number of EOFS: 10,5,5 for the 3 input variables
#' data <- prepareData(x = x, y = y, spatial.predictors = list(n.eofs = c(10,5,5)))
#' # Using joined PCs as predictors. Explained variance 95%
#' data <- prepareData(x = x, y = y, 
#' spatial.predictors = list(v.exp = 0.95, which.combine =getVarNames(x)))
#' # Using local predictors: the 4 closest gridboxes
#' data <- prepareData(x = x, y = y,local.predictors = list(n=4, vars = getVarNames(x)))
#' # Using joined PCs and local predictors: the 4 closest gridboxes
#' data <- prepareData(x = x, y = y,local.predictors = list(n=4, vars = getVarNames(x)),
#' spatial.predictors = list(v.exp = 0.95, which.combine =getVarNames(x)))

prepareData <- function(x, y, global.vars = NULL, combined.only = TRUE, spatial.predictors = NULL, local.predictors = NULL, extended.predictors = NULL) {
    y <- getTemporalIntersection(obs = y, prd = x, which.return = "obs")
    x <- getTemporalIntersection(obs = y, prd = x, which.return = "prd")
    dates <- getRefDates(x)
    xyCoords <- getCoordinates(x)
    
    # LOCAL PREDICTORS
    nn <- NULL
    all.predictor.vars <- if (!is.null(global.vars)) {
        global.vars
    } else {
        getVarNames(x)
    }
    if (!is.null(local.predictors)) {
        all.predictor.vars <- c(all.predictor.vars, local.predictors$vars)
        pred.nn.ind <- predictor.nn.indices(vars = local.predictors$vars,
                                            n = local.predictors$n,
                                            x = x, y = y) 
        local.predictors[["local.index.list"]] <- pred.nn.ind
        nn <- predictor.nn.values(nn.indices.list = pred.nn.ind,
                                  grid = x,
                                  fun = local.predictors$fun)
        attr(nn, "local.index.list") <- pred.nn.ind
    }
    all.predictor.vars %<>% unique()
    # GLOBAL PREDICTOR SUBSETTING
    if (!is.null(global.vars)) {
        x  %<>%  subsetGrid(var = global.vars, drop = FALSE)
    }
    varnames <- getVarNames(x)
    # spatial predictors
    full.pca <- NULL
    if (!is.null(spatial.predictors)) {
        spatial.predictors[["grid"]] <- x
        if (is.null(spatial.predictors$which.combine)) {
            message("NOTE: The COMBINED PC won't be calculated as 'which.combine' argument in spatial.predictors is missing (with no default)")
        }
        full.pca <- do.call("prinComp", spatial.predictors)
        x <- if (!is.null(spatial.predictors$which.combine)) {
            if (length(varnames) == length(spatial.predictors$which.combine)) combined.only <- TRUE
            if (combined.only) {
                  if (!is.null(full.pca$COMBINED)) {
                        full.pca$COMBINED[[1]]$PCs
                  }else if (length(varnames) == 1) {
                        full.pca[[varnames]][[1]]$PCs
                  }
            } else {
                combined.vars <- attributes(full.pca$COMBINED)$combined_variables
                all.vars <- names(full.pca)[-length(full.pca)]
                aux <- full.pca[all.vars[which(!all.vars %in% combined.vars)]]
                aux <- lapply(1:length(aux), function(i) aux[[i]][[1]]$PCs)
                cbind(do.call("cbind", aux), full.pca$COMBINED[[1]]$PCs)
            }
        } else {
            aux <- lapply(1:length(full.pca), function(i) full.pca[[i]][[1]]$PCs)
            do.call("cbind", aux)
        }
    # RAW predictors    
    } 
    else {
        # Construct a predictor matrix of raw input grids
        x <- redim(x, var = TRUE)
        nvars <- getShape(x, "var")
        suppressMessages(
            aux.list <- lapply(1:nvars, function(i) {
                  subsetGrid(x, var = varnames[i]) %>% redim(drop = TRUE) %>% extract2("Data") %>% array3Dto2Dmat()
            })
        )
        x <- do.call("cbind", aux.list)
        aux.list <- NULL
        
        if (!is.null(extended.predictors$n)) {
          inp.n <- ncol(x) + 1
          hid.n <- extended.predictors$n
          r <- 2*((inp.n) ** (-0.5))
          w <- array(data = runif(inp.n*hid.n, min = -r, max = r), c(inp.n,hid.n))
          if (!is.null(extended.predictors$module)) {
            ww <- w[1:(inp.n - 1),]
            dim(ww) <- c(7,7,20,hid.n)
            mask <- array(data = 0,dim = dim(ww))
            for (zzz in 1:hid.n) {
              r1 <- sample(1:(dim(ww)[1] - extended.predictors$module),size = 1)
              r2 <- sample(1:(dim(ww)[2] - extended.predictors$module),size = 1)
              mask[(r1:(r1 + extended.predictors$module)),(r2:(r2 + extended.predictors$module)),,zzz] <- ww[(r1:(r1 + extended.predictors$module)),(r2:(r2 + extended.predictors$module)),,zzz]
            }
            dim(mask) <- c(980,hid.n)
            w[1:(inp.n - 1),] <- mask
          }
          w[inp.n,] <- runif(hid.n, min = -2, max = 2)
          x.bias <- cbind(x,array(data = 1,dim = c(nrow(x), 1))) 
          h <- x.bias %*% w
          x <- 1/(1 + exp(-h))
        }
    }
    predictor.list <- list("y" = y, "x.global" = x, "x.local" = nn, "pca" = full.pca)
    if (!is.null(spatial.predictors)) spatial.predictors <- spatial.predictors[-grep("grid", names(spatial.predictors))]
    if (!is.null(extended.predictors$n)) attr(predictor.list,"auxRandomMatrix") <- w
    if (!is.null(spatial.predictors) & is.null(local.predictors)) {
      nature <- "spatial"
    } else if (is.null(spatial.predictors) & !is.null(local.predictors)) {
      nature <- "local"
    } else if (!is.null(spatial.predictors) & !is.null(local.predictors)) {
      nature <- "spatial+local"
    } else if (!is.null(extended.predictors)) {
      nature <- "extended"
    } else if (is.null(spatial.predictors) & is.null(local.predictors) & is.null(extended.predictors)) {
      nature <- "raw"
    }
    attr(predictor.list, "spatialPred.pars") <- spatial.predictors
    attr(predictor.list, "localPred.pars") <- local.predictors
    attr(predictor.list, "predictor.vars") <- all.predictor.vars
    attr(predictor.list, "nature") <- nature
    attr(predictor.list, "globalPred.vars") <- varnames
    attr(predictor.list, "dates") <- dates
    attr(predictor.list, "xyCoords") <- xyCoords
    
    return(predictor.list)
}


#' @title Calculate spatial index position of nearest neighbors
#' @description Calculate spatial index position of nearest neighbors
#' @param vars Character vector of the names of the variable(s) to be used as local predictors
#' @param n Integer vector. Number of nearest neighbours to be considered for each variables in \code{vars}.
#' Its length must be either 1 --in case there is one single local predictor variable or the same number of neighbours is 
#' desired for every local predictor variable-- or equal the length of \code{vars}, to consider a varying number of 
#' neighbours for each local predictor variable.
#' @param x.grid The (multi)grid containing the predictors
#' @param y.grid The grid containing the predictand
#' @keywords internal
#' @importFrom transformeR getVarNames 
#' @return A list of matrices (rows = local neighbors, cols = predictand points). The list is of the same
#' length as the number of local predictor variables chosen
#' @family downscaling.helpers
#' @author J Bedia

predictor.nn.indices <- function(vars = NULL, n = NULL, x, y) {
      if (is.null(vars)) stop("Undefined local predictor variables. A value for 'vars' argument is required", call. = FALSE)
      if (is.null(n)) stop("Undefined number of local neighbours. A value for 'n' argument is required", call. = FALSE)
      varnames <- getVarNames(x)
      if (!any(vars %in% varnames)) stop("The requested neigh.var was not found in the predictor grid", call. = FALSE)
      # Selection of nearest neighbour indices
      # Coordinates matrices
      coords.y <- get2DmatCoordinates(y)
      coords.x <- get2DmatCoordinates(x)
      if (any(n > nrow(coords.x))) stop("Too many neighbours selected (more than predictor grid cells)", call. = FALSE) 
      if (length(n) == 1) {
            n <- rep(n, length(vars))
      }
      if (length(n) != length(vars)) {
            stop("Incorrect number of neighbours selected: this should be either 1 or a vector of the same length as 'vars'", call. = FALSE)
      }
      # The index list has the same length as the number of local predictor variables, containing a matrix of index positions
      local.pred.list <- lapply(1:length(vars), function(j) {
            ind.mat <- vapply(1:nrow(coords.y), FUN.VALUE = numeric(n[j]), FUN = function(i) {
                  dists <- sqrt((coords.y[i,1] - coords.x[,1])^2 + (coords.y[i,2] - coords.x[,2])^2)
                  which(dists %in% sort(dists)[1:n[j]])
            })
      })
      names(local.pred.list) <- vars
      return(local.pred.list)
}

#' @title Construct local predictor matrices
#' @description Constructs the local predictor matrices given their spatial index position 
#' (as returned by \code{\link{predictor.nn.indices}}) and additional parameters.
#' @param nn.indices.list A list of index positions as returned by \code{\link{predictor.nn.indices}}
#' @param fun A list of funtions for neighbour aggregation for each local predictor variable.
#'  Default to \code{NULL}, and no aggregation is performed
#' @param grid The grid containing the local predictor information. It can be either the predictors (training)
#'  grid or the simulation (prediction) grid, depending on the situation.
#' @return A nested list of 2D matrices with the following structure: sites/members
#' @keywords internal
#' @family downscaling.helpers
#' @importFrom transformeR subsetGrid array3Dto2Dmat
#' @importFrom magrittr %>% extract2
#' @author J Bedia

predictor.nn.values <- function(nn.indices.list, grid, fun = NULL) {
    # Aggregation function of neighbors
    if (length(nn.indices.list) == 1) {
        fun <- list(fun)
    } 
    if (!is.null(fun) && (length(fun) != length(nn.indices.list))) {
        stop("Incorrect number of neighbours selected: this should be either 1 or a vector of the same length as 'nn.indices.list'", call. = FALSE)
    }
    out.list <- lapply(1:length(nn.indices.list), function(j) {
        vars <- names(nn.indices.list)
        aux <- subsetGrid(grid, var = vars[j], drop = TRUE) %>% redim(member = TRUE) 
        n.mem <- getShape(aux, "member")
        mem.list <- lapply(1:n.mem, function(i) {
            aux1 <- subsetGrid(aux, members = i, drop = TRUE) %>% extract2("Data") %>% array3Dto2Dmat()
            # Local predictor matrix list
            l <- lapply(1:ncol(nn.indices.list[[j]]), function(k) aux1[ ,nn.indices.list[[j]][ ,k]])
            # Neighbour aggregation function 
            if (!is.null(fun[[j]])) {
                arg.list <- c(fun[[j]], "MARGIN" = 1)
                lapply(1:length(l), function(k) do.call(what = "apply", c("X" = l[k], arg.list))) 
            } else {
                l
            }
        })
        names(mem.list) <- paste("member", 1:n.mem, sep = "_")
        return(mem.list)
    })
    names(out.list) <- names(nn.indices.list)
    lapply(1:length(out.list[[1]][[1]]), function(i) {
        aux <- lapply(1:length(out.list[[1]]), function(j) {
            expr <- paste0("cbind(", paste0("out.list[[", 1:length(out.list), "]][[", j,"]][[",i,"]]", collapse = ","), ")")
            parse(text = expr) %>% eval()
        })
        names(aux) <- paste("member", 1:length(aux), sep = "_")
        return(aux)
    })
 }


