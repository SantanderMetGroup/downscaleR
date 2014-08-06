#if (length(dimDiff)==0){
#  F <- calibrateProj(obs$Data, aperm(pred$Data,dimPerm), aperm(sim$Data,dimPerm), method = method, varcode = obs$Variable$varName, pr.threshold = threshold)
#  sim$Data<-aperm(F,c(dimPerm[2:length(dim(obs$Data))],dimPerm[1]))
#}else{
#  P<-aperm(pred$Data,dimPerm)
#  F<-aperm(sim$Data,dimPerm)
#}
#if (length(dimDiff)==1){
#  for (n in 1:dim(sim$Data)[length(dim(obs$Data))+1]){
#    F1 <- calibrateProj(obs$Data, P[,,,n], F[,,,n], method = method, varcode = obs$Variable$varName, pr.threshold = threshold)
#    sim$Data[,,,n]<-aperm(F1,c(dimPerm[2:length(dim(obs$Data))],dimPerm[1]))
#  }
#}
#if (length(dimDiff)==2){
#  for (n in 1:dim(sim$Data)[length(dim(obs$Data))+1]){
#    for (m in 1:dim(sim$Data)[length(dim(obs$Data))+2]){
#      F1 <- calibrateProj(obs$Data, P[,,,n,m], F[,,,n,m], method = method, varcode = obs$Variable$varName, pr.threshold = threshold)
#      sim$Data[,,,n,m]<-aperm(F1,c(dimPerm[2:length(dim(obs$Data))],dimPerm[1]))
#    }
#  }
#}
#if (length(dimDiff)==3){
#  for (n in 1:dim(sim$Data)[length(dim(obs$Data))+1]){
#    for (m in 1:dim(sim$Data)[length(dim(obs$Data))+2]){
#      for (t in 1:dim(sim$Data)[length(dim(obs$Data))+3]){
#        F1 <- calibrateProj(obs$Data, P[,,,n,m,t], F[,,,n,m,t], method = method, varcode = obs$Variable$varName, pr.threshold = threshold)
#        sim$Data[,,,n,m,t]<-aperm(F1,c(dimPerm[2:length(dim(obs$Data))],dimPerm[1]))
#      }
#    }
#  }
#}
#if (any(grepl(obs$Variable$varName,c("pr","tp","precipitation","precip")))){
#  attr(sim$Data, "threshold") <-  threshold
#}
#attr(sim$Data, "correction") <-  method
#return(sim)
#}

#################################################################################################################
# ISI-MIP Bias correction method (http://www.pik-potsdam.de/research/climate-impacts-and-vulnerabilities/research/rd2-cross-cutting-activities/isi-mip/about/isi-mip-fast-track)
# Citation: Hempel, S., Frieler, K., Warszawski, L., Schewe, J., and Piontek, F.: A trend-preserving bias correction ? the ISI-MIP approach, Earth Syst. Dynam., 4, 219-236, doi:10.5194/esd-4-219-2013, 2013.
#----------------------------------------------------------------------------------
#' @title Bias correction methods
#' @description Implementation of several standard bias correction methods
#' @author S. Herrera \email{sixto@@predictia.es}
#' @export
#' 
isimip <- function (obs, pred, sim, pr.threshold = 1) {

# threshold<-pr.threshold
# if size(threshold(:),1)==1,threshold=[1 1 1]*threshold;end

datesObs <- as.POSIXct(obs$Dates$start, tz="GMT", format="%Y-%m-%d %H:%M:%S")
datesObs<-cut(datesObs, "month")
yearList<-unlist(strsplit(as.character(datesObs), "[-]"))
months<-unique(unlist(strsplit(as.character(datesObs), "[-]"))[seq(2,length(yearList),3)])
dimObs<-dim(obs$Data)
dimObs[which(grepl("^time",attr(obs$Data, "dimensions")))]<-length(unique(datesObs))
dimPred<-dim(pred$Data)
dimPred[which(grepl("^time",attr(pred$Data, "dimensions")))]<-length(unique(datesObs))
monthlyObs<- array(data = NA, dim = dimObs)
monthlyPred<- array(data = NA, dim = dimPred)
month2dayObs<- array(data = NA, dim = c(length(datesObs),1))
obs.time.index <- grep("^time$", attr(obs$Data, "dimensions"))
pred.time.index <- grep("^time$", attr(pred$Data, "dimensions"))
for (i in 1:length(unique(datesObs))){
  indMonth<-which(datesObs == unique(datesObs)[i])
  month2dayObs[indMonth]<-i
  indTimeObs <- rep(list(bquote()), length(dimObs))
  for (d in 1:length(dimObs)){
    indTimeObs[[d]] <- 1:dimObs[d]
  }
  indTimeObs[[obs.time.index]] <- i
  indTimeObs<-as.matrix(expand.grid(indTimeObs))
  indTimeObs1 <- rep(list(bquote()), length(dimObs))
  indTimeObs1[[obs.time.index]] <- indMonth
  callObs <- as.call(c(list(as.name("["),quote(obs$Data)), indTimeObs1))
  monthlyObs[indTimeObs] <- apply(eval(callObs), FUN = mean, MARGIN = setdiff(1:length(dimObs),obs.time.index), na.rm = TRUE)
  indTimeObs <- rep(list(bquote()), length(dimPred))
  for (d in 1:length(dimPred)){
    indTimeObs[[d]] <- 1:dimPred[d]
  }
  indTimeObs[[pred.time.index]] <- i
  indTimeObs<-as.matrix(expand.grid(indTimeObs))
  indTimeObs1 <- rep(list(bquote()), length(dimPred))
  indTimeObs1[[pred.time.index]] <- indMonth
  callObs <- as.call(c(list(as.name("["),quote(pred$Data)), indTimeObs1))
  monthlyPred[indTimeObs] <- apply(eval(callObs), FUN = mean, MARGIN = setdiff(1:length(dimPred),pred.time.index), na.rm = TRUE)
}
datesFor <- as.POSIXct(sim$Dates$start, tz="GMT", format="%Y-%m-%d %H:%M:%S")
datesFor<-cut(datesFor, "month")
dimFor<-dim(sim$Data)
dimFor[which(grepl("^time",attr(sim$Data, "dimensions")))]<-length(unique(datesFor))
monthlyFor<- array(data = NA, dim = dimFor)
month2dayFor<- array(data = NA, dim = c(length(datesFor),1))
sim.time.index <- grep("^time$", attr(sim$Data, "dimensions"))
for (i in 1:length(unique(datesFor))){
  indMonth<-which(datesFor == unique(datesFor)[i])
  month2dayFor[indMonth]<-i
  indTimeObs <- rep(list(bquote()), length(dimFor))
  for (d in 1:length(dimFor)){
    indTimeObs[[d]] <- 1:dimFor[d]
  }
  indTimeObs[[sim.time.index]] <- i
  indTimeObs<-as.matrix(expand.grid(indTimeObs))
  indTimeObs1 <- rep(list(bquote()), length(dimFor))
  indTimeObs1[[sim.time.index]] <- indMonth
  callObs <- as.call(c(list(as.name("["),quote(sim$Data)), indTimeObs1))
  monthlyFor[indTimeObs] <- apply(eval(callObs), FUN = mean, MARGIN = setdiff(1:length(dimFor),sim.time.index), na.rm = TRUE)
}
if (any(grepl(obs$Variable$varName,c("tas","mean temperature")))){
  dimAux<-dimPred
  dimAux[pred.time.index]<-length(months)
  monthlyCorrection<-array(data = NA, dim = dimAux)
  for (i in 1:length(months)){
    indMonth<-which(grepl(paste("-",months[i],"-", sep=""),unique(datesObs)))
    indTimeObs <- rep(list(bquote()), length(dimObs))
    indTimeObs[[obs.time.index]] <- indMonth
    callObs <- as.call(c(list(as.name("["),quote(monthlyObs)), indTimeObs))
    auxObs <- apply(eval(callObs), FUN = mean, MARGIN = setdiff(1:length(dimObs),obs.time.index), na.rm = TRUE)
    indTimeObs <- rep(list(bquote()), length(dimPred))
    for (d in 1:length(dimPred)){
      indTimeObs[[d]] <- 1:dimPred[d]
    }
    indTimeObs<-as.matrix(expand.grid(indTimeObs))
    indTimeObs1 <- rep(list(bquote()), length(dimPred))
    indTimeObs1[[pred.time.index]] <- indMonth
    callObs <- as.call(c(list(as.name("["),quote(monthlyPred)), indTimeObs1))
    auxPrd <- apply(eval(callObs), FUN = mean, MARGIN = setdiff(1:length(dimPred),pred.time.index), na.rm = TRUE)
    indTimeObs <- rep(list(bquote()), length(dimPred))
    for (d in 1:length(dimPred)){
      indTimeObs[[d]] <- 1:dimPred[d]
    }
    indTimeObs[[pred.time.index]] <- i
    indTimeObs<-as.matrix(expand.grid(indTimeObs))
    monthlyCorrection[indTimeObs]<-array(auxObs, dim = dim(auxPrd))-auxPrd
  }
# Replicar lo de los indices para cada una de las opciones  
  
  pred$Data <- pred$Data-monthlyPred[month2dayObs,,,]
  obs$Data <- obs$Data-monthlyObs[,,month2dayObs]
  sim$Data <- sim$Data-monthlyFor[month2dayFor,,,]
  for (i in 1:dim(obs$Data)[1]){
    for (j in 1:dim(obs$Data)[2]){
      for (m in 1:length(months)){
        indMonth<-which(grepl(paste("-",months[m],"-", sep=""),unique(datesObs)))
        indMonthFor<-which(grepl(paste("-",months[m],"-", sep=""),unique(datesFor)))
        if (any(!is.na(obs$Data[i,j,indMonth]))){
          for (k in 1:dim(pred$Data)[4]){
            if (any(!is.na(pred$Data[indMonth,i,j,k]))){
              lmAdjust<-lm(sort(pred$Data[indMonth,i,j,k], decreasing = FALSE, na.last = NA) ~ sort(obs$Data[i,j,indMonth], decreasing = FALSE, na.last = NA)-1)
              sim$Data[indMonthFor,i,j,k]<-coef(lmAdjust)[1]*sim$Data[indMonthFor,i,j,k]
            }
          }
        }
      }
    }
  }
  sim$Data <- sim$Data+monthlyFor[month2dayFor,,,]
  for (i in 1:length(months)){
    indMonthFor<-which(grepl(paste("-",months[i],"-", sep=""),unique(datesFor)))
    if (!is.null(indMonthFor)){
      for (j in 1:length(indMonthFor)){
        sim$Data[indMonthFor[j],,,] <- sim$Data[indMonthFor[j],,,]+monthlyCorrection[i,,,]
      }
    }
  }
}
if (any(grepl(obs$Variable$varName,c("pr","tp","precipitation","precip")))){
  attr(sim$Data, "threshold") <-  threshold
}
attr(sim$Data, "correction") <-  "ISI-MIP"

return(sim)
}

################################################################################################
#dimDiff<-NULL
#dimPerm <- 1:length(attr(pred$Data, "dimensions"))
#if (length(setdiff(attr(pred$Data, "dimensions"),attr(obs$Data, "dimensions")))>0){
#  dimDiff<-setdiff(attr(pred$Data, "dimensions"),attr(obs$Data, "dimensions"))
#  for (k in 1:length(dimDiff)) {
#    dimPerm[which(grepl(dimDiff[k],attr(pred$Data, "dimensions")))]<-length(attr(obs$Data, "dimensions"))+k
#  }
#}
#dimPerm[which(dimPerm<=length(dim(obs$Data)))]<-  match(attr(obs$Data, "dimensions"), attr(pred$Data, "dimensions"))
#
#isimip <- function (obs, pred, sim, datesObs = 1:dim(obs)[1] , datesFor = 1:dim(sim)[1] , origin = NULL,
#    tg = NULL, wss = NULL, threshold = 0, variable = c("pr","rss","rsds","rls","rlds","ps","wss","huss","hus","tas","tasmax","tasmin","uas","vas","ua","va")) {
#    obs <- obs
#    pred <- pred
#    sim <- sim
    # En este paso se eliminan instancias con valores perdidos de cada una de las series (4 marzo)
#    c("obs","pred","sim") -> arg.names
#    c() -> ind
#    for (i in 1:length(arg.names)) {
#        vec <- get(arg.names[i])
#        if (any(is.na(vec))) {
#            a <- which(is.na(vec))
#            ind <- c(ind, a)
#        }
#    }
    #print(ind)
#    if (length(ind) > 0) {
#        obs <- obs[-ind]
#        pred <- pred[-ind]
#        sim <- sim[-ind]
#    }
    # Fin cambio
#	ndata<-dim(obs)[1]
#	Nest<-dim(obs)[2]
	# Agregacion mensual:
#	datesObs<-chron(datesObs, format = "yyyy-mm-dd", origin = origin)
#	daysCut<-seq(datesObs[1], datesObs[ndata], by = "months")
#	monthlyO<-matrix(data = NA, nrow = length(daysCut)[1], ncol = Nest)
#	monthlyP<-matrix(data = NA, nrow = length(daysCut)[1], ncol = Nest)
#	k1<-length(daysCut)[1]-1
#	for (k in 1:k1) {
#		ind1<-which(datesObs==daysCut[k])
#		ind2<-which(datesObs==daysCut[k+1])-1
#		monthlyO[k,] <- mean(obs[ind1:ind2,], na.rm = TRUE)
#		monthlyP[k,] <- mean(pred[ind1:ind2,], na.rm = TRUE)
#	}
#	ind1<-which(datesObs==daysCut[length(daysCut)[1]])
#	ind2<-ndata
#	monthlyO[length(daysCut)[1],] <- mean(obs[ind1:ind2,], na.rm = TRUE)
#	monthlyP[length(daysCut)[1],] <- mean(pred[ind1:ind2,], na.rm = TRUE)
#	meses=unique(months(datesObs))

	
#	nbins<-dim(sim)[1]
#    obs <- sort(obs)
#    pred <- sort(pred)
#    ix <- sort(sim, index.return = TRUE)$ix
#    sim <- sort(sim)
#	auxQ<-1:nbins
#	auxQ<-100*auxQ/nbins
#	auxQ<-seq(1/nbins, 1, by=1/nbins) 
#	if (ndata==nbins){
#		delta<-sim-pred
#    delta_i <- sim - pred
#    delta <- delta_i - mean(delta_i)
#	}
#	else{
#		aux<-quantile(pred, auxQ, na.rm = TRUE, names = FALSE)
#		delta<-sim-aux
#	}
#	deltaMean<-mean(sim, na.rm = TRUE)-mean(pred, na.rm = TRUE)
#    varcode <- match.arg(varcode, c("tas", "hurs", "pr", "wss"))
#	if (varcode == "pr") {
#		nP<-matrix(0,1,Nest)
#		Pth<-matrix(0,1,Nest)
#		for (k in 1:Nest) {
#			nP[k]<-length(which(obs[,k]<=threshold))
#			if (pred[nP[k]+1,k]<=threshold){
#				ind<-which(pred[,k]>=threshold)
#				params<-params<-fitdistr(obs[nP(k)+1:ind[1],k],"gamma"),# shape<-params$estimate[1]; scale<-1/params$estimate[2]
#				pred(nP(k)+1:ind[1],k)<-rgamma(length(nP(k)+1:ind[1]),params$estimate[1],params$estimate[2]);
#				pred[,k] <- sort(pred[,k])
#			}
#			pred[1:min(nP[k],ndata),k] <- 0
#			Pth[k] <- pred[nP[k]+1,k]
#		}
#	}
#    method <- match.arg(method, c("qqadj", "qqmap", "bias"))
#	normfun <- match.arg(normfun, c("prctile", "std"))
#    return.par <- return.par
#    prj <- rep(NA, nbins, Nest)
#    if (method == "bias") {
#		aux<-quantile(obs, auxQ, na.rm = TRUE, names = FALSE)
#		for (k in 1:Nest) {
#			prj[,k]<-aux[,k]+deltaMean[k]
#		}
#    }
#    if (method == "biasratio") {
#		for (k in 1:Nest) {
#			prj[,k]<-sim[,k]*(mean(obs[,k],na.rm=TRUE)/mean(prd[,k],na.rm=TRUE))
#		}
#    }
#    if (method == "biaslog") {
#		for (k in 1:Nest) {
#			oo<-obs[,k]
#			oo[which(obs[,k]<1e-16)]<-0
#			pp<-pred[,k]
#			pp[which(pred[,k]<1e-16)]<-0
#			prj[,k]<-sim[,k]^(mean(log(obs[,k]),na.rm=TRUE)/mean(log(prd[,k]),na.rm=TRUE))
#		}
#    }
#    if (method == "qqmap") {
#		aux<-quantile(obs, auxQ, na.rm = TRUE, names = FALSE)
#		for (k in 1:Nest) {
#			prj[,k]<-aux[,k]+delta[,k]
#		}
#    }
#    if (method == "qqadj") {
#        if (varcode == "tas") {
#			g <- matrix(1,1,Nest)
#        }
#        else {
#            g <- mean(obs, na.rm = TRUE)/mean(pred, na.rm = TRUE)
#        }
#        if (varcode == "pr") {
#            obs[which(obs>threshold)] -> obsn
#            pred[which(pred>threshold)] -> predn
#			if (normfun == "prctile"){
#				f <- (quantile(obsn, 0.9, names = FALSE) - quantile(obsn, 
#					0.1, names = FALSE))/(quantile(predn, 0.9, names = FALSE) - 
#					quantile(predn, 0.1, names = FALSE))
#			} else{
#				f <- sd(obsn, na.rm = TRUE)/sd(predn, na.rm = TRUE)
#			}
#            ########
#        }
#        else {
#			if (normfun == "prctile"){
#				f <- (quantile(obs, 0.75, names = FALSE) - quantile(obs, 
#					0.25, names = FALSE))/(quantile(pred, 0.75, names = FALSE) - 
#					quantile(pred, 0.25, names = FALSE))
#			} else{
#				f <- sd(obsn, na.rm = TRUE)/sd(predn, na.rm = TRUE)
#			}
#        }
#    }
#    prj <- rep(NA, length(obs))
#    for (k in 1:length(obs)) {
#        prj[k] <- obs[k] + g * mean(delta_i) + f * delta[k]
#    }
#    if (varcode == "pr") {
#        ############## Modificacion (29Feb2012)
#        if (length(which(pred == 0)) > 0) {
#            nzp <- length(which(sim == 0)) * length(which(obs == 0))/length(which(pred == 0))
#            if (nzp > length(obs) | is.finite(nzp) == FALSE) {
#                length(obs) -> nzp
#            }
#            prj[1:nzp] <- 0
#        }
#        ############# Fin cambio
#    }
#    prix <- cbind(ix, prj)
#    prj <- prix[order(prix[, 1]), 2]
#    if (return.par == TRUE) {
#        prj <- list(corrvals = prj, g = g, f = f, Mean_delta = mean(delta_i))
#    }
#    return(prj)
#}
# End
#----------------------------------------------------------------------------------