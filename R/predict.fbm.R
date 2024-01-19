#' Title
#'
#' @import stats
#'
#' @keywords internal
predict.fbm <- function(omega,time.app,time.test, h){
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
      pred[i] <- omega[which(time.app==time_inf)]*0.5*(time.test[i]^(2*h)+max(time.app)^(2*h)-abs(time.test[i]-max(time.app))^(2*h))/(max(time.app)^(2*h))
      }
    }
    if (length(sup)>0 & length(inf)==0){
      time_sup <- min(time.app[sup])
      pred[i] <- omega[which(time.app==time_sup)]*0.5*(time.test[i]^(2*h)+min(time.app)^(2*h)-abs(time.test[i]-min(time.app))^(2*h))/(min(time.app)^(2*h))
    }
    return(pred)
  }
}
