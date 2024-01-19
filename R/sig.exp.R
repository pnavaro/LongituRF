#' Title
#'
#' @param sigma
#' @param id
#' @param Z
#' @param epsilon
#' @param Btilde
#' @param time
#' @param sigma2
#' @param alpha
#'
#' @import stats
#'
#' @keywords internal
sig.exp <- function(Y,sigma,id,Z, epsilon, Btilde, time, sigma2,alpha){ #### fonction d'actualisation du param?tre de la variance des erreurs
  nind <- length(unique(id))
  Nombre <- length(id)
  sigm <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    K <- cov.exp(time[w], alpha)
    V <- Z[w,]%*%Btilde%*%t(Z[w,])+diag(rep(sigma,length(Y[w])))+sigma2*K
    sigm <- sigm + t(epsilon[w])%*%epsilon[w] + sigma*(length(w)-sigma*(sum(diag(solve(V)))))
  }
  sigm <- sigm/Nombre
  return(sigm)
}
