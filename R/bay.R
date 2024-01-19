#' Title
#'
#' @import stats
#'
#' @keywords internal
bay <- function(bhat,Bhat,Z,id, sigmahat){ #### actualisation des param?tres de B
  nind <- length(unique(id))
  q <- dim(Z)[2]
  Nombre <- length(id)
  D <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    V <- Z[w,, drop=FALSE]%*%Bhat%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigmahat),length(w),length(w))
    D <- D+ (bhat[j,]%*%t(bhat[j,]))+ (Bhat- Bhat%*%t(Z[w,, drop=FALSE])%*%solve(V)%*%Z[w,, drop=FALSE]%*%Bhat)
  }
  D <- D/nind
  return(D)
}
