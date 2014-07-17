#' @title Bias correction methods
#' @description Implementation of several standard bias correction methods
#' @author S. Herrera \email{sixto@@predictia.es}
#' @export


biasCorrection <- function (obs, pred, sim, method = c("qqmap", "delta", "unbiasing", "piani"), varcode = c("tas", "hurs", "pr", "wss"), pr.threshold = 1) {
	obs <- obs
	pred <- pred
	sim <- sim
	varcode <- match.arg(varcode, c("tas","hurs","pr","wss"))
	if (varcode == "pr") {
		threshold<-pr.threshold
		nP<-matrix(data = NA, ncol=dim(pred)[3], nrow=dim(pred)[2])
		Pth<-matrix(data = NA, ncol=dim(pred)[3], nrow=dim(pred)[2])
		for (i in 1:dim(pred)[2]){
			for (j in 1:dim(pred)[3]){
				nP[i,j]<-sum(as.double(obs[,i,j]<=threshold & !is.na(obs[,i,j])), na.rm = TRUE)
				if (nP[i,j]>0){
					ix<-sort(pred[,i,j], decreasing = FALSE, na.last = NA, index.return = TRUE)$ix
					Ps<-sort(pred[,i,j], decreasing = FALSE, na.last = NA)
					Pth[i,j]<-Ps[nP[i,j]+1]
					if (Ps[nP[i,j]+1]<=threshold){
						ind<-min(which(Ps > threshold & !is.na(Ps)))
						# [Shape parameter  Scale parameter]
						Os<-sort(obs[,i,j], decreasing = FALSE, na.last = NA)
						auxGamma<-fitdistr(Os[(nP[i,j]+1):ind],"gamma")
						Ps[(nP[i,j]+1):ind]<-rgamma(ind-nP[i,j], auxGamma$estimate[1], rate = auxGamma$estimate[2])
						Ps<-sort(Ps, decreasing = FALSE, na.last = NA)
					}
					ind<-min(nP[i,j],dim(pred)[1])
					Ps[1:ind]<-0
					pred[ix,i,j]<-Ps
				}
			}
		}
	}
	method <- match.arg(method, c("qqmap","scaling","unbiasing","piani"))
	if (method == "unbiasing") {
		obsMean <- apply(obs, FUN = mean, MARGIN = c(2,3), na.rm = TRUE)
		prdMean <- apply(pred, FUN = mean, MARGIN = c(2,3), na.rm = TRUE)
		for (k in 1:dim(sim)[1]) {
			sim[k,,]<-sim[k,,]-prdMean+obsMean
		}
	}
	if (method == "scaling") {
		obsMean <- apply(obs, FUN = mean, MARGIN = c(2,3), na.rm = TRUE)
		prdMean <- apply(pred, FUN = mean, MARGIN = c(2,3), na.rm = TRUE)
		for (k in 1:dim(sim)[1]) {
			sim[k,,]<-(sim[k,,]/prdMean)*obsMean
		}
	}
	if (method == "qqmap") {
		if (varcode != "pr") {	
			for (i in 1:dim(sim)[2]) {
				for (j in 1:dim(sim)[3]) {
					if (any(!is.na(pred[,i,j])) & any(!is.na(obs[,i,j]))) {
						ePrd<-ecdf(pred[,i,j])
						sim[,i,j]<-quantile(obs[,i,j], probs = ePrd(sim[,i,j]), na.rm = TRUE, type = 8)
					}
				}
			}
		} else {
			for (i in 1:dim(sim)[2]) {
				for (j in 1:dim(sim)[3]) {
					if (any(!is.na(pred[,i,j])) & any(!is.na(obs[,i,j]))){
						ePrd<-ecdf(pred[which(pred[,i,j]>Pth[i,j]),i,j])
						noRain<-which(sim[,i,j]<=Pth[i,j] & !is.na(sim[,i,j]))
						rain<-which(sim[,i,j]>Pth[i,j] & !is.na(sim[,i,j]))
						drizzle<-which(sim[,i,j]>Pth[i,j] & sim[,i,j]<=min(pred[which(pred[,i,j]>Pth[i,j]),i,j], na.rm = TRUE) & !is.na(sim[,i,j]))
						eFrc<-ecdf(sim[rain,i,j])
						sim[drizzle,i,j]<-quantile(sim[which(sim[,i,j]>min(pred[which(pred[,i,j]>Pth[i,j]),i,j], na.rm = TRUE) & !is.na(sim[,i,j])),i,j], probs = eFrc(sim[drizzle,i,j]), na.rm = TRUE, type = 8)
						sim[rain,i,j]<-quantile(obs[which(obs[,i,j]>threshold & !is.na(obs[,i,j])),i,j], probs = ePrd(sim[rain,i,j]), na.rm = TRUE, type = 8)
						sim[noRain,i,j]<-0
					}
				}
			}
		}
	}
	if (method == "piani" & varcode == "pr") {
		for (i in 1:dim(sim)[2]){
			for (j in 1:dim(sim)[3]){
				if (nP[i,j]>0){
					ind<-which(obs[,i,j]>threshold & !is.na(obs[,i,j]))
					obsGamma<-fitdistr(obs[ind,i,j],"gamma")
					ind<-which(pred[,i,j]>0 & !is.na(pred[,i,j]))
					prdGamma<-fitdistr(pred[ind,i,j],"gamma")
					rain<-which(sim[,i,j]>Pth[i,j] & !is.na(sim[,i,j]))
					noRain<-which(sim[,i,j]<=Pth[i,j] & !is.na(sim[,i,j]))
					auxF<-pgamma(sim[rain,i,j],prdGamma$estimate[1], rate = prdGamma$estimate[2])
					sim[rain,i,j]<-qgamma(auxF,obsGamma$estimate[1], rate = obsGamma$estimate[2])
					sim[noRain,i,j]<-0
				}
			}
		}
	}
	return(sim)
}
