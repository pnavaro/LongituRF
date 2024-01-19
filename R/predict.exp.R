#' blabla
#'
#' @import stats
#'
#' @keywords internal
predict.exp <- function(omega,time.app,time.test, alpha){
  pred <- rep(0,length(time.test))
  for (i in 1:length(time.test)){
    inf <- which(time.app<=time.test[i])
    sup <- which(time.app>time.test[i])
    if (length(inf)>0){
      if (length(sup)>0){
        time_inf <- max(time.app[inf])
        time_sup <- min(time.app[sup])
        pred[i] <- mean(c(omega[which(time.app==time_inf)],omega[which(time.app==time_sup)]))}
      if(length(sup)==0) {time_inf <- max(time.app[inf])
      pred[i] <- omega[which(time.app==time_inf)]*exp(-alpha*abs(time.test[i]-max(time.app)))
      }
    }
    if (length(sup)>0 & length(inf)==0){
      time_sup <- min(time.app[sup])
      pred[i] <- omega[which(time.app==time_sup)]*exp(-alpha*abs(time.test[i]-max(time.app)))
    }
    return(pred)
  }
}
