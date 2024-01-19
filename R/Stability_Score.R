#' Stability score function for (S)MERF and (S)REEMforest methods
#'
#' Computes the stability scores for (S)MERF and (S)REEMforest methods.
#'
#' @param X [matrix]: A \code{N}x\code{p} matrix containing the \code{p} predictors of the fixed effects, column codes for a predictor.
#' @param Y [vector]: A vector containing the output trajectories.
#' @param id [vector]: Is the vector of the identifiers for the different trajectories.
#' @param Z [matrix]: A \code{N}x\code{q} matrix containing the \code{q} predictor of the random effects.
#' @param time [vector]: Is the vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' @param mtry [numeric]: Number of variables ramdomly picked to split each node.
#' @param ntree [numeric]: Number of trees in the RF.
#' @param sto [string]: Defines the covariance function of the stochastic process, can be either \code{"none"} for no stochastic process, \code{"BM"} for Brownian motion, \code{OrnUhl} for standard Ornstein-Uhlenbeck process, \code{BBridge} for Brownian Bridge, \code{fbm} for Fractional Brownian motion; can also be a function defined by the user.
#' @param method [string]: Defines the method to be used, can be either "MERF" or "REEMforest".
#' @param eta [numeric]: The size of the neighborhood for the stability score. Can be a vector, in this case, returns the stability scores corresponding to all the values of the vector.
#' @param nvars [numeric]: The number of variables to consider among the most impotant variables. Can be a vector, in this case, the function returns the stability scores corresponding to all the values of the vector.
#'
#' @export
#'
#' @return A matrix with all the stability scores corresponding to the eta and nvars values. The $i$th row corresponds to the $i$th value of eta while the $i$th column corresponds to the $i$ value of nvars.
#
Stability_Score <- function(X,Y,Z,id,time,mtry,ntree, sto="BM",method="MERF", eta = c(1:ncol(X)),nvars=c(1:ncol(X))){

  if (method=="REEMforest"){
    sortie1 <- REEMforest(X=X,Y=Y,Z=Z,id=id,time=time,mtry=mtry,ntree=ntree,sto=sto)
    sortie2 <- REEMforest(X=X,Y=Y,Z=Z,id=id,time=time,mtry=mtry,ntree=ntree,sto=sto)
  }

  else {
    sortie1 <- MERF(X=X,Y=Y,Z=Z,id=id,time=time,mtry=mtry,ntree=ntree,sto=sto)
    sortie2 <- MERF(X=X,Y=Y,Z=Z,id=id,time=time,mtry=mtry,ntree=ntree,sto=sto)
  }

  imp1 <- sort(sortie1$forest$importance[,1], decreasing = TRUE, index.return=TRUE)
  imp2 <- sort(sortie2$forest$importance[,1], decreasing = TRUE, index.return=TRUE)

  ss <- matrix(NA,length(eta),length(nvars))
  for (i in 1:length(eta)){
    for (k in 1:length(nvars)){
      nind = rep(0,nvars[k])
      for (l in 1:nvars[k]){
        variable = imp1$ix[l]
        variable2_interv = imp2$ix[c(max(1,l-eta[i]):min(nvars[k],l+eta[i]))]
        nind[l] = 1*length(which(variable2_interv == variable)>0)
      }
      ss[i,k] = mean(nind)
    }
  }
  SS <- as.data.frame(ss)
  colnames(SS) = nvars
  rownames(SS) = eta
  return(SS)
}
