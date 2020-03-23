#   prepareNewData.R Configuration of data for downscaling method predictions
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
#' @param data.structure A structure, as returned by \code{\link{prepareData}}
#' @return A named list with the components required by the downscaling method in order to perform the predictions
#' @export
#' @seealso \href{https://github.com/SantanderMetGroup/downscaleR/wiki/preparing-predictor-data}{downscaleR Wiki} for preparing predictors for downscaling and seasonal forecasting.
#' @author J Bedia
#' @family downscaling.helpers
#' @importFrom transformeR getVarNames subsetGrid redim getShape getCoordinates grid2PCs getRefDates array3Dto2Dmat grid2PCs
#' @importFrom magrittr %>% extract2 
#' @examples
#' # Loading data
#' data("VALUE_Iberia_tas")
#' y <- VALUE_Iberia_tas 
#' data("NCEP_Iberia_hus850", "NCEP_Iberia_psl", "NCEP_Iberia_ta850")
#' x <- makeMultiGrid(NCEP_Iberia_hus850, NCEP_Iberia_psl, NCEP_Iberia_ta850)
#' # Calculating EOFs
#' data <- prepareData(x = x, y = y, spatial.predictors = list(v.exp = 0.95))
#' # Projecting a new dataset to the calculated EOFs
#' newdata <- prepareNewData(x,data)

prepareNewData <- function(newdata, data.structure) {
  x.varnames <- attr(data.structure, "predictor.vars")
  newdata.varnames <- getVarNames(newdata)
  if (!all(x.varnames %in% newdata.varnames)) stop("The variable names in predictor and newdata structures do not match", call. = FALSE)
  # Spatial check primero para segurar que los nn.indices son correctos
  if (!identical(attr(data.structure, "xyCoords"), getCoordinates(newdata))) stop("Spatial mismatch between predictor and newdata", call. = FALSE)
  # Local predictors 
  newdata.local.list <- NULL # A list containing, for each location of the predictand, the local neighbour values
  if (!is.null(data.structure$x.local)) {
    local.index.list <- attributes(data.structure$x.local)$"local.index.list"
    newdata.local.list <- predictor.nn.values(nn.indices.list = local.index.list,
                                              grid = newdata)
    attr(newdata.local.list,"predictorNames") <- attributes(data.structure$x.local)$"predictorNames"
  }
  # Global predictors
  global.pred.vars <- attr(data.structure, "globalPred.vars")
  ## Raw
  if (is.null(data.structure$pca)) {
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
    if (!is.null(attr(data.structure,"auxRandomMatrix"))) {
      newdata.global.list <- lapply(1:length(newdata.global.list), function(z){
        h <- cbind(newdata.global.list[[z]],rep(1,nrow(newdata.global.list[[z]]))) %*% attr(data.structure,"auxRandomMatrix")
        1/(1 + exp(-h))
      })
    }
    names(newdata.global.list) <- paste("member", 1:length(newdata.global.list), sep = "_")
    attr(newdata.global.list,"predictorNames") <- attributes(data.structure$x.global)$"predictorNames"
  } else {
    ## PCA predictors   
    if (is.null(data.structure$pca$COMBINED)) { # The COMBINED PC is not being used
      pca.mat.list <- lapply(1:length(global.pred.vars), function(i) {
        aux <- subsetGrid(newdata, var = global.pred.vars[i]) %>% grid2PCs(prinCompObj = data.structure$pca)
      })
      names(pca.mat.list) <- global.pred.vars
      newdata.global.list <- lapply(1:length(pca.mat.list[[1]]), function(x) {
        expr <- paste0("cbind(", paste0("pca.mat.list[[", 1:length(pca.mat.list), "]][[x]]", collapse = ","), ")")
        parse(text = expr) %>% eval()
      })
      pca.mat.list <- NULL
    } else {# The COMBINED PC is being used
      combined.pred.vars <- attr(data.structure$pca$COMBINED, "combined_variables")
      pca.mat.list <- lapply(1:length(combined.pred.vars), function(i) {
        z <- match(combined.pred.vars[i],global.pred.vars)
        aux <- subsetGrid(newdata, var = global.pred.vars[z]) %>% redim(member = TRUE)
        escala <- attributes(data.structure$pca[[z]][[1]]$orig)$'scaled:scale'
        center <- attributes(data.structure$pca[[z]][[1]]$orig)$'scaled:center'
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
        combined.std.mat[[i]] %*% data.structure$pca$COMBINED[[1]]$EOFs
      })
      combined.std.mat <- NULL
    }
    names(newdata.global.list) <- paste("member", 1:length(newdata.global.list), sep = "_")
    attr(newdata.global.list,"predictorNames") <- attributes(data.structure$x.global)$"predictorNames"
  }
  newdata.refdates <- list(start = getRefDates(newdata, "start"), end = getRefDates(newdata, "end"))
  attr(newdata.global.list,"nature") <- attr(data.structure,"nature")
  return(list("x.global" = newdata.global.list,
              "x.local" = newdata.local.list,
              "Dates" = newdata.refdates)
  )
}

