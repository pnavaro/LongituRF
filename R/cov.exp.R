#' Title
#'
#' @param time
#' @param alpha
#'
#' @import stats
#'
#' @keywords internal
cov.exp <- function(time,alpha){
  K <- matrix(0,length(time),length(time))
  for (i in 1:length(time)){
    for (j in 1:length(time)){
      K[i,j] <- exp(-(alpha*abs(time[i]-time[j])))
    }
  }
  return(K)
}
