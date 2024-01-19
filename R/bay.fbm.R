#' Title
#'
#' @param bhat
#' @param Bhat
#' @param Z
#' @param id
#' @param sigmahat
#' @param time
#' @param sigma2
#' @param h
#'
#' @import stats
#'
#' @keywords internal
bay.fbm <- function(bhat,Bhat,Z,id, sigmahat, time, sigma2,h){ #### actualisation des param?tres de B
  nind <- length(unique(id))
  q <- dim(Z)[2]
  Nombre <- length(id)
  D <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    K <- cov.fbm(time[w], h)
    V <- Z[w,]%*%Bhat%*%t(Z[w,])+diag(rep(sigmahat,length(w)))+sigma2*K
    D <- D+ (bhat[j,]%*%t(bhat[j,]))+ (Bhat- Bhat%*%t(Z[w,])%*%solve(V)%*%Z[w,]%*%Bhat)
  }
  D <- D/nind
  return(D)
}
