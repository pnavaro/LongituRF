#' Title
#'
#'
#' @import stats
#'
#' @keywords internal
Moy_fbm <- function(id,Btilde,sigmahat,Phi,Y,Z, H, time, sigma2){
  S1<- 0
  S2<- 0
  nind <- length(unique(id))
  for (i in 1:nind){
    w <- which(id==unique(id)[i])
    K <- cov.fbm(time[w],H)
    V <- Z[w,, drop=FALSE]%*%Btilde%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigmahat),length(w),length(w))+ sigma2*K
    S1 <- S1 + t(Phi[w,, drop=FALSE])%*%solve(V)%*%Phi[w,, drop=FALSE]
    S2 <- S2 + t(Phi[w,, drop=FALSE])%*%solve(V)%*%Y[w]
  }
  M <- solve(S1)%*%S2
}
