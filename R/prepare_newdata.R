#   prepare_newdata.R Configuration of data for downscaling method predictions
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

#' @title Prepare newdata for predictions
#' @description Prepare the prediction data according to the definition of the experiment
#' @param newdata A grid containing the prediction data.
#' @param predictor A predictor structure, as returned by \code{\link{prepare_predictors}}
#' @return A named list with the components required by the downscaling method in order to perform the predictions
#' @export
#' @author J Bedia
#' @examples 
#' # See the dedicated vignette by typing:
#' # utils::vignette(topic = "configuring_newdata", package = "downscaleR")
#' @family downscaling.helpers
#' @importFrom transformeR getVarNames subsetGrid redim getShape getCoordinates grid2PCs getRefDates array3Dto2Dmat grid2PCs
#' @importFrom magrittr %>% extract2 

prepare_newdata <- function(newdata, predictor) {
    x.varnames <- attr(predictor, "predictor.vars")
    newdata.varnames <- getVarNames(newdata)
    if (!all(x.varnames %in% newdata.varnames)) stop("The variable names in predictor and newdata structures do not match", call. = FALSE)
    # Spatial check primero para segurar que los nn.indices son correctos
    if (!identical(attr(predictor, "xyCoords"), getCoordinates(newdata))) stop("Spatial mismatch between predictor and newdata", call. = FALSE)
    # Local predictors 
    newdata.local.list <- NULL # A list containing, for each location of the predictand, the local neighbour values
    if (!is.null(predictor$x.local)) {
        local.index.list <- attributes(predictor$x.local)$"local.index.list"
        newdata.local.list <- predictor.nn.values(nn.indices.list = local.index.list,
                                                  grid = newdata)
    }
    # Global predictors
    global.pred.vars <- attr(predictor, "globalPred.vars")
    ## Raw
    if (is.null(predictor$pca)) {
        # A list (must handle members)
        newdata.global.list <- lapply(1:length(global.pred.vars), function(i) {
            aux <- subsetGrid(newdata, var = global.pred.vars[i]) %>% redim(var = FALSE, member = TRUE) 
            n.mem <- getShape(aux, "member")
            l <- lapply(1:n.mem, function(j) {
                subsetGrid(aux, members = j, drop = TRUE) %>% extract2("Data") %>% array3Dto2Dmat()
            })
            names(l) <- paste("member", 1:n.mem, sep = "_")
            return(l)
        })
        names(newdata.global.list) <- global.pred.vars
        newdata.global.list <- lapply(1:length(newdata.global.list[[1]]), function(i) {
            expr <- paste0("cbind(", paste0("newdata.global.list[[", 1:length(newdata.global.list), "]][[i]]", collapse = ","), ")")
            parse(text = expr) %>% eval()
        })
        if (!is.null(attr(predictor,"neurons.matrix"))) {
          newdata.global.list <- lapply(1:length(newdata.global.list), function(z){
            h <- cbind(newdata.global.list[[z]],rep(1,nrow(newdata.global.list[[z]]))) %*% attr(predictor,"neurons.matrix")
            1/(1 + exp(-h))
          })
        }
        names(newdata.global.list) <- paste("member", 1:length(newdata.global.list), sep = "_")
        
    } else {
    ## PCA predictors   
        if (is.null(predictor$pca$COMBINED)) { # The COMBINED PC is not being used
            pca.mat.list <- lapply(1:length(global.pred.vars), function(i) {
                aux <- subsetGrid(newdata, var = global.pred.vars[i]) %>% grid2PCs(prinCompObj = predictor$pca)
            })
            names(pca.mat.list) <- global.pred.vars
            newdata.global.list <- lapply(1:length(pca.mat.list[[1]]), function(x) {
                expr <- paste0("cbind(", paste0("pca.mat.list[[", 1:length(pca.mat.list), "]][[x]]", collapse = ","), ")")
                parse(text = expr) %>% eval()
            })
            pca.mat.list <- NULL
        } else {# The COMBINED PC is being used
            pca.mat.list <- lapply(1:length(global.pred.vars), function(i) {
                aux <- subsetGrid(newdata, var = global.pred.vars[i]) %>% redim(member = TRUE)
                escala <- attributes(predictor$pca[[i]][[1]]$orig)$'scaled:scale'
                center <- attributes(predictor$pca[[i]][[1]]$orig)$'scaled:center'
                n.mem <- getShape(aux, "member")
                lapply(1:n.mem, function(j) {
                    subsetGrid(aux, members = j, drop = TRUE) %>% extract2("Data") %>% array3Dto2Dmat() %>% scale(center = center, scale = escala)
                })
            })
            ## The combined standardized matrix is created, as a list of members: 
            combined.std.mat <- lapply(1:length(pca.mat.list[[1]]), function(i) {
                expr <- paste0("cbind(", paste0("pca.mat.list[[", 1:length(pca.mat.list), "]][[", i,"]]", collapse = ","), ")")
                parse(text = expr) %>% eval()
            })
            pca.mat.list <- NULL
            # The COMBINED PCs are calculated using the predictor EOFs:
            newdata.global.list <- lapply(1:length(combined.std.mat), function(i) {
                combined.std.mat[[i]] %*% predictor$pca$COMBINED[[1]]$EOFs
            })
            combined.std.mat <- NULL
        }
        names(newdata.global.list) <- paste("member", 1:length(newdata.global.list), sep = "_")
    }
    newdata.refdates <- list(start = getRefDates(newdata, "start"), end = getRefDates(newdata, "end"))
    return(list("x.global" = newdata.global.list,
                "x.local" = newdata.local.list,
                "Dates" = newdata.refdates)
    )
}

