#' Predict with longitudinal trees and random forests.
#'
#' @param object : a \code{longituRF} output of (S)MERF; (S)REEMforest; (S)MERT or (S)REEMtree function.
#' @param X [matrix]: matrix of the fixed effects for the new observations to be predicted.
#' @param Z [matrix]: matrix of the random effects for the new observations to be predicted.
#' @param id [vector]: vector of the identifiers of the new observations to be predicted.
#' @param time [vector]: vector of the time measurements of the new observations to be predicted.
#' @param ... : low levels arguments.
#'
#' @import stats
#' @import randomForest
#'
#' @return vector of the predicted output for the new observations.
#'
#' @export
#'
#' @examples \donttest{
#' set.seed(123)
#' data <- DataLongGenerator(n=20) # Generate the data composed by n=20 individuals.
#' REEMF <- REEMforest(X=data$X,Y=data$Y,Z=data$Z,id=data$id,time=data$time,mtry=2,ntree=500,sto="BM")
#' # Then we predict on the learning sample :
#' pred.REEMF <- predict(REEMF, X=data$X,Z=data$Z,id=data$id, time=data$time)
#' # Let's have a look at the predictions
#' # the predictions are in red while the real output trajectories are in blue:
#' par(mfrow=c(4,5),mar=c(2,2,2,2))
#' for (i in unique(data$id)){
#'   w <- which(data$id==i)
#'   plot(data$time[w],data$Y[w],type="l",col="blue")
#'   lines(data$time[w],pred.REEMF[w], col="red")
#' }
#' # Train error :
#' mean((pred.REEMF-data$Y)^2)
#'
#' # The same function can be used with a fitted SMERF model:
#' smerf <-MERF(X=data$X,Y=data$Y,Z=data$Z,id=data$id,time=data$time,mtry=2,ntree=500,sto="BM")
#' pred.smerf <- predict(smerf, X=data$X,Z=data$Z,id=data$id, time=data$time)
#' # Train error :
#' mean((pred.smerf-data$Y)^2)
#' # This function can be used even on a MERF model (when no stochastic process is specified)
#' merf <-MERF(X=data$X,Y=data$Y,Z=data$Z,id=data$id,time=data$time,mtry=2,ntree=500,sto="none")
#' pred.merf <- predict(merf, X=data$X,Z=data$Z,id=data$id, time=data$time)
#' # Train error :
#' mean((pred.merf-data$Y)^2)
#'
#' }
predict.longituRF <- function(object, X,Z,id,time,...){
  n <- length(unique(id))
  id_btilde <- object$id_btilde
  f <- predict(object$forest,X)
  Time <- object$time
  id_btilde <- object$id_btilde
  Ypred <- rep(0,length(id))
  id.app=object$id
  if (object$sto=="none"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      k <- which(id_btilde==unique(id)[i])
      Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,]
    }
    return(Ypred)
  }

  if (object$sto=="exp"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      k <- which(id_btilde==unique(id)[i])
      om <- which(id.app==unique(id)[i])
      Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,] + predict.exp(object$omega[om],Time[om],time[w], object$alpha)
    }
    return(Ypred)
  }

  if (object$sto=="fbm"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      k <- which(id_btilde==unique(id)[i])
      om <- which(id.app==unique(id)[i])
      Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,] + predict.fbm(object$omega[om],Time[om],time[w], object$Hurst)
    }
    return(Ypred)
  }

  for (i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    k <- which(id_btilde==unique(id)[i])
    om <- which(id.app==unique(id)[i])
    Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,] + predict.sto(object$omega[om],Time[om],time[w], object$sto)
  }
  return(Ypred)
}
