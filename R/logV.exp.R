#' Title
#'
#' @param Y
#' @param Z
#' @param id
#' @param fhat
#' @param sigmahat
#' @param sigma2
#' @param H
#' @param B
#'
#' @import stats
#'
#' @keywords internal
logV.exp <- function(Y, fhat, Z, time, id, B, sigma2,sigmahat,H){
  logl <- 0
  for (i in 1:length(unique(id))){
    indiv <- which(id==unique(id)[i])
    K <- cov.exp(time[indiv], H)
    V <- Z[indiv,]%*%B%*%t(Z[indiv,])+ sigma2*K+ as.numeric(sigmahat)*diag(rep(1,length(indiv)))
    logl <- logl + log(det(V)) + t(Y[indiv]-fhat[indiv])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
  }
  return(-logl)
}
