#' Title
#'
#' @param time
#' @param H
#'
#' @keywords internal
cov.fbm <- function(time,H){
  K <- matrix(0,length(time),length(time))
  for (i in 1:length(time)){
    for (j in 1:length(time)){
      K[i,j] <- 0.5*(time[i]^(2*H)+time[j]^(2*H)-abs(time[i]-time[j])^(2*H))
    }
  }
  return(K)
}
