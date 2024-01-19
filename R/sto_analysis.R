#' Title
#'
#' @import stats
#'
#' @keywords internal
sto_analysis <- function(sto, time){
  MAT <- matrix(0,length(time), length(time))

  if (class(sto)=="function"){
    for (i in 1:length(time)){
      for (j in 1:length(time)){
        MAT[i,j] <- sto(time[i], time[j])
      }
    }
    return(MAT)
  }

  if (sto=="BM"){
    for (i in 1:length(time)){
      for (j in 1:length(time)){
        MAT[i,j] <- min(time[i], time[j])
      }
    }
    return(MAT)
  }

  if (sto=="OrnUhl"){
    for (i in 1:length(time)){
      for (j in 1:length(time)){
        MAT[i,j] <- exp(-abs(time[i]-time[j])/2)
      }
    }
    return(MAT)
  }

  if (sto=="BBridge"){
    for (i in 1:length(time)){
      for (j in 1:length(time)){
        MAT[i,j] <- min(time[i], time[j]) - time[i]*time[j]
      }
    }
    return(MAT)
  }

}
