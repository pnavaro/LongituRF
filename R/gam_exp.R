#' Title
#'
#' @param sigm
#' @param id
#' @param Z
#' @param B
#' @param time
#' @param sigma2
#' @param omega
#' @param id_omega
#' @param alpha
#'
#' @import stats
#'
#' @keywords internal
gam_exp <- function(Y,sigm,id,Z,Btilde,time,sigma2,omega,id_omega, alpha){
  nind <- length(unique(id))
  Nombre <- length(id)
  gam <- 0
  for (k in 1:nind){
    indiv <- which(id==unique(id)[k])
    K <- cov.exp(time[indiv],alpha)
    V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,,drop=FALSE])+diag(rep(sigm,length(Y[indiv])))+ sigma2*K
    Omeg <- omega[k,which(id_omega[k,]==1)]
    gam <-gam+ (t(Omeg)%*%solve(K)%*%Omeg) + sigma2*(length(indiv)-sigma2*sum(diag(solve(V)%*%K)))
  }
  return(as.numeric(gam)/Nombre)
}
