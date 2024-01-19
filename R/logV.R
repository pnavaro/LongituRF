#' Title
#'
#'
#' @import stats
#'
#' @keywords internal
logV <- function(Y,f,Z,time,id,B,gamma,sigma, sto){
  Vraisem <- 0
  if (sto=="none"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      V <- Z[w,,drop=FALSE]%*%B%*%t(Z[w,,drop=FALSE])+diag(as.numeric(sigma),length(w),length(w))
      print(V)
      print(det(V))
      V2 = log(det(V))+ t(Y[w]-f[w])%*%solve(V)%*%(Y[w]-f[w])
      print(V2)
      if (V2<Inf){
        Vraisem <- Vraisem + V2
      }
    }
    return(Vraisem)
  }
  for (i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    K <- sto_analysis(sto,time[w])
    V <- Z[w,,drop=FALSE]%*%B%*%t(Z[w,,drop=FALSE])+gamma*K+ diag(as.numeric(sigma),length(w),length(w))
    Vraisem <- Vraisem + log(det(V))+ t(Y[w]-f[w])%*%solve(V)%*%(Y[w]-f[w])
  }
  return(Vraisem)
}
