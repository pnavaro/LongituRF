#' Title
#'
#' @import stats
#'
#' @keywords internal
predict.sto <- function(omega,time.app,time.test, sto){
  pred <- rep(0,length(time.test))

  if (class(sto)=="function"){
    for (i in 1:length(time.test)){
      inf <- which(time.app<=time.test[i])
      sup <- which(time.app>time.test[i])
      if (length(inf)>0){
        if (length(sup)>0){
          time_inf <- max(time.app[inf])
          time_sup <- min(time.app[sup])
          pred[i] <- mean(c(omega[which(time.app==time_inf)],omega[which(time.app==time_sup)]))}
        if(length(sup)==0) {time_inf <- max(time.app[inf])
        pred[i] <- omega[which(time.app==time_inf)]*((sto(time.test[i],max(time.app)))/sto(max(time.app),max(time.app)))
        }
      }
      if (length(sup)>0 & length(inf)==0){
        time_sup <- min(time.app[sup])
        pred[i] <- omega[which(time.app==time_sup)]*((sto(time.test[i],min(time.app)))/sto(min(time.app),min(time.app)))
      }
    }
    return(pred)
  }

  else {
    for (i in 1:length(time.test)){
      inf <- which(time.app<=time.test[i])
      sup <- which(time.app>time.test[i])
      if (length(inf)>0){
        if (length(sup)>0){
          time_inf <- max(time.app[inf])
          time_sup <- min(time.app[sup])
          pred[i] <- mean(c(omega[which(time.app==time_inf)],omega[which(time.app==time_sup)]))}
        if(length(sup)==0) {time_inf <- max(time.app[inf])
        if (sto=="BM"){
          pred[i] <- omega[which(time.app==time_inf)]}
        if (sto=="OrnUhl"){
          pred[i] <- omega[which(time.app==time_inf)]*(exp(-abs(time.test[i]-max(time.app))/2))
        }
        if (sto=="BBridge"){
          pred[i] <- omega[which(time.app==time_inf)]*((1-time.test[i])/(1-max(time.app)^2))
        }
        }
      }
      if (length(sup)>0 & length(inf)==0){
        time_sup <- min(time.app[sup])
        if (sto=="BM"){
          pred[i] <- omega[which(time.app==time_sup)]*(time.test[i]/min(time.app))}
        if (sto=="OrnUhl"){
          pred[i] <- omega[which(time.app==time_sup)]*(exp(-abs(time.test[i]-min(time.app))/2))
        }
        if (sto=="BBridge"){
          pred[i] <- omega[which(time.app==time_sup)]*(time.test[i]/min(time.app))
        }
      }
    }}
  return(pred)
}
