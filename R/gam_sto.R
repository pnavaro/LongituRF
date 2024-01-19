#' Title
#'
#' @import stats
#'
#' @keywords internal
gam_sto <- function(sigma,id,Z, Btilde, time, sigma2,sto, omega){
  nind <- length(unique(id))
  Nombre <- length(id)
  gam <- 0
  for (k in 1:nind){
    indiv <- which(id==unique(id)[k])
    K <- sto_analysis(sto,time[indiv])
    V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,,drop=FALSE])+diag(as.numeric(sigma),length(indiv),length(indiv))+ sigma2*K
    Omeg <- omega[indiv]
    gam <-gam+ (t(Omeg)%*%solve(K)%*%Omeg) + sigma2*(length(indiv)-sigma2*sum(diag(solve(V)%*%K)))
  }
  return(as.numeric(gam)/Nombre)
}
