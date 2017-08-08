# load("ignore/juaco/data/obsPredSim_NCEP.Rdata", verbose = TRUE)
# # obs <- obs.tmean
# 
# data("VALUE_Iberia_tas")
# obs <- VALUE_Iberia_tas
# 
# getShape(pred)
# getShape(sim)
# 
# # Check the consistency of x and newdata
# checkDim(pred, sim, dimensions =  c("var","lon", "lat"))
# 
# 
# # Obtain the coordinates of x (=newdata) and y
# coords.y <- get2DmatCoordinates(obs)
# coords.x <- get2DmatCoordinates(pred)
# 
# # Compute PCA
# xpc <- prinComp(pred, v.exp = .9)
# 
# # Target variable(s) of neighbors (either number or name)
# neigh.vars = c("psl", "tas")
# n.neigh = c(9,3) # number of neighbors
# neigh.fun = list(list(FUN = "mean", na.rm = TRUE), NULL)

#' @param neigh.vars Character vectorwith the names of the variable(s ) to be used as local predictors
#' @param n.neighs Integer vector. Number of nearest neighbours to be considered for each variables in \code{neigh.vars}.
#' Its length must be either 1 --in case there is one single local predictor variable or the same number of neighbours is 
#' desired for every local predictor variable-- or equal the length of \code{neigh.vars}, to consider a varying number of 
#' neighbours for each local predictor variable.
#' @param x.grid The (multi)grid containing the predictors
#' @param y.grid The grid containing the predictand
#' @keywords internal
#' @importFrom transformeR getVarNames subsetGrid array3Dto2Dmat
#' @importFrom magrittr %>% extract2
#' @return A list of matrices (rows = local neighbors, cols = predictand points). The ist is of the same
#' length as the number of local predictor variables chosen
#' @family downscaling.helpers
#' @author J Bedia

predictor.nn.indices <- function(neigh.vars, n.neighs, x.grid, y.grid) {
    varnames <- getVarNames(pred)
    varnames <- c("hus850", varnames)
    if (!any(neigh.vars %in% varnames)) stop("The requested neigh.var was not found in the predictor grid", call. = FALSE)
    # Selection of nearest neighbour indices
    if (any(n.neigh > nrow(coords.x))) stop("Too many neighbours selected (more than predictor grid cells)", call. = FALSE) 
    if (length(n.neigh) == 1) {
        n.neigh <- rep(n.neigh, length(neigh.vars))
    } else if (length(n.neigh) != length(neigh.vars)) {
        stop("Incorrect number of neighbours selected: this should be either 1 or a vector of the same length as 'neigh.vars'", call. = FALSE)
    }
    # Coordinates matrices
    coords.y <- get2DmatCoordinates(y.grid)
    coords.x <- get2DmatCoordinates(x.grid)
    # The index list has the same length as the number of local predictor variables, containing a matrix of index positions
    local.pred.list <- lapply(1:length(neigh.vars), function(j) {
        ind.mat <- vapply(1:nrow(coords.y), numeric(n.neigh[j]), FUN = function(i) {
            dists <- sqrt((coords.y[i,1] - coords.x[,1])^2 + (coords.y[i,2] - coords.x[,2])^2)
            which(dists %in% sort(dists)[1:n.neigh[j]])
        })
    })
    names(local.pred.list) <- neigh.vars
    return(local.pred.list)
}

a <- predictor.nn.indices(neigh.vars = c("psl", "tas"), n.neighs = c(10,3), x.grid = pred, y.grid = obs)



#' @param nn.indices.list
#' @param neigh.fun A list of funtions for neighbour aggregation for each local predictor variable.
#'  Default to \code{NULL}, and no aggregation is performed
#' @param grid The grid containing the local predictor information. It can be either the predictors (training)
#'  grid or the simulation (prediction) grid, depending on the situation.
#' @keywords internal
#' @family downscaling.helpers
#' @author J Bedia

#  nn.indices.list = a
#  str(a)
# grid = pred
# neigh.fun = list(NULL, list(FUN = mean, na.rm = TRUE))

predictor.nn.values <- function(nn.indices.list, neigh.fun, grid) {
    # Aggregation function of neighbors
    if (length(nn.indices.list) == 1) {
        neigh.fun <- list(neigh.fun)
    } else if (length(neigh.fun) != length(nn.indices.list)) {
        stop("Incorrect number of neighbours selected: this should be either 1 or a vector of the same length as 'nn.indices.list'", call. = FALSE)
    }
    out.list <- lapply(1:length(nn.indices.list), function(j) {
        aux <- subsetGrid(grid, var = neigh.vars[j], drop = TRUE) %>% extract2("Data") %>% array3Dto2Dmat()
        # Local predictor matrix list
        l <- lapply(1:ncol(nn.indices.list[[j]]), function(k) aux[,nn.indices.list[[j]][,k]])
        # Neighbour aggregation function 
        if (!is.null(neigh.fun[[j]])) {
            arg.list <- c(neigh.fun[[j]], "MARGIN" = 1)
            lapply(1:length(l), function(k) do.call(what = "apply", c("X" = l[k], arg.list))) # %>% as.matrix(ncol = 1))
        } else {
            l
        }
    })
    lapply(1:length(out.list[[1]]), function(x) cbind(out.list[[1]][[x]], out.list[[2]][[x]]))
}

# q <- predictor.nn.values(a, neigh.fun = list(NULL, list(FUN = mean, na.rm = TRUE)), grid = pred)