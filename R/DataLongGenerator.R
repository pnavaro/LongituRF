#' Longitudinal data generator
#'
#'
#' Simulate longitudinal data according to the semi-parametric stochastic mixed-effects model given by: \deqn{Y_i(t)=f(X_i(t))+Z_i(t)\beta_i + \omega_i(t)+\epsilon_i}
#' with \eqn{Y_i(t)} the output at time \eqn{t} for the \eqn{i}th individual; \eqn{X_i(t)} the input predictors (fixed effects) at time \eqn{t} for the \eqn{i}th individual;
#' \eqn{Z_i(t)} are the random effects at time \eqn{t} for the \eqn{i}th individual; \eqn{\omega_i(t)} is a Brownian motion with volatility \eqn{\gamma^2=0.8} at time \eqn{t} for the \eqn{i}th individual; \eqn{\epsilon_i} is the residual error with
#' variance \eqn{\sigma^2=0.5}.
#' The data are simulated according to the simulations in low dimensional in the low dimensional scheme of the paper <doi:10.1177/0962280220946080>
#'
#' @param n [numeric]: Number of individuals. The default value is \code{n=50}.
#' @param p [numeric]: Number of predictors. The default value is \code{p=6}.
#' @param G [numeric]: Number of groups of predictors with temporal behavior, generates \code{p-G} input variables with no temporal behavior.
#'
#' @import mvtnorm
#' @import latex2exp
#'
#' @return a list of the following elements: \itemize{
#' \item \code{Y:} vector of the output trajectories.
#' \item \code{X :} matrix of the fixed-effects predictors.
#' \item \code{Z:} matrix of the random-effects predictors.
#' \item \code{id: } vector of the identifiers for each individual.
#' \item \code{time: } vector the the time measurements for each individual.
#' }
#'
#' @export
#'
#' @examples
#' oldpar <- par()
#' oldopt <- options()
#' data <- DataLongGenerator(n=17, p=6,G=6) # Generate the data
#' # Let's see the output :
#' w <- which(data$id==1)
#' plot(data$time[w],data$Y[w],type="l",ylim=c(min(data$Y),max(data$Y)), col="grey")
#' for (i in unique(data$id)){
#'   w <- which(data$id==i)
#'   lines(data$time[w],data$Y[w], col='grey')
#' }
#' # Let's see the fixed effects predictors:
#' par(mfrow=c(2,3), mar=c(2,3,3,2))
#' for (i in 1:ncol(data$X)){
#'   w <- which(data$id==1)
#'   plot(data$time[w],data$X[w,i], col="grey",ylim=c(min(data$X[,i]),
#'   max(data$X[,i])),xlim=c(1,max(data$time)),main=latex2exp::TeX(paste0("$X^{(",i,")}$")))
#'   for (k in unique(data$id)){
#'     w <- which(data$id==k)
#'     lines(data$time[w],data$X[w,i], col="grey")
#'   }
#' }
#' par(oldpar)
#' options(oldopt)
#'
DataLongGenerator <- function(n=50,p=6,G=6){

  mes <-floor(4*runif(n)+8)
  time <- NULL
  id <- NULL
  nb2 <- c(1:n)
  for (i in 1:n){
    time <- c(time, seq(1,mes[i], by=1))
    id <- c(id, rep(nb2[i], length(seq(1,mes[i], by=1))))
  }

  bruit <- floor(0*p)
  bruit <- bruit+ (p-bruit)%%G
  nices <- NULL
  for (i in 1:G){
    nices <- c(nices,rep(i,(p-bruit)/G))
  }

  comportements <- matrix(0,length(time),G)
  comportements[,1] <- 2.44+0.04*(time-((time-6)^2)/(time/3))
  comportements[,2] <- 0.5*time-0.1*(time-5)^2
  comportements[,3] <- 0.25*time-0.05*(time-6)^2
  comportements[,4] <- cos((time-1)/3)
  comportements[,5] <- 0.1*time + sin(0.6*time+1.3)
  comportements[,6] <- -0.1*time^2


  X <- matrix(0,length(time), p)
  for (i in 1:(p-bruit)){
    X[,i] <- comportements[,nices[i]] + rnorm(length(time),0 ,0.2)
  }

  for (j in 1:n){
    w <- which(id==j)
    X[w,1:(p-bruit)] <- X[w,1:(p-bruit)] + rnorm(1,0,0.1)
  }

  for (i in (p-bruit):p){
    X[,i] <- rnorm(length(time),0, 3)
  }

  f <- 1.3*X[,1]^2 + 2*sqrt(abs(X[,which(nices==2)[1]]))

  sigma <- cbind(c(0.5,0.6),c(0.6,3))
  Btilde<- matrix(0,length(unique(id)),2)
  for (i in 1:length(unique(id))){
    Btilde[i,] <- rmvnorm(1, mean=rep(0,2),sigma=sigma)
  }

  Z <- as.matrix(cbind(rep(1,length(f)),2*runif(length(f))))

  effets  <- NULL
  for (i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    effets <- c(effets, Z[w,, drop=FALSE]%*%Btilde[i,])
  }
  ##### simulation de mouvemments brownien
  gam <- 0.8
  BM <- NULL
  m <- length(unique(id))
  for (i in 1:m){
    w <- which(id==unique(id)[i])
    W <- rep(0,length(w))
    t <- time[w]
    for (j in 2:length(w)){
      W[j] <- W[j-1]+sqrt(gam*(t[j]-t[j-1]))*rnorm(1,0,1)
    }
    BM <- c(BM,W)
  }

  sigma2 <- 0.5
  Y <- f + effets +rnorm(length(f),0,sigma2)+BM
  return(list(Y=Y,X=X,Z=Z,id=id, time=time))
}
