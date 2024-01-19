#' Title
#'
#' @import stats
#'
#' @keywords internal
Moy <- function(id,Btilde,sigmahat,Phi,Y,Z){
  S1<- 0
  S2<- 0
  nind <- length(unique(id))
  for (i in 1:nind){
    w <- which(id==unique(id)[i])
    V <- Z[w,, drop=FALSE]%*%Btilde%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigmahat),length(w),length(w))
    S1 <- S1 + t(Phi[w,, drop=FALSE])%*%solve(V)%*%Phi[w,, drop=FALSE]
    S2 <- S2 + t(Phi[w,, drop = FALSE])%*%solve(V)%*%Y[w]
  }
  return(solve(S1)%*%S2)
}
