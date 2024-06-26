% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataLongGenerator.R
\name{DataLongGenerator}
\alias{DataLongGenerator}
\title{Longitudinal data generator}
\usage{
DataLongGenerator(n = 50, p = 6, G = 6)
}
\arguments{
\item{n}{\link{numeric}: Number of individuals. The default value is \code{n=50}.}

\item{p}{\link{numeric}: Number of predictors. The default value is \code{p=6}.}

\item{G}{\link{numeric}: Number of groups of predictors with temporal behavior, generates \code{p-G} input variables with no temporal behavior.}
}
\value{
a list of the following elements: \itemize{
\item \code{Y:} vector of the output trajectories.
\item \code{X :} matrix of the fixed-effects predictors.
\item \code{Z:} matrix of the random-effects predictors.
\item \code{id: } vector of the identifiers for each individual.
\item \code{time: } vector the the time measurements for each individual.
}
}
\description{
Simulate longitudinal data according to the semi-parametric stochastic mixed-effects model given by: \deqn{Y_i(t)=f(X_i(t))+Z_i(t)\beta_i + \omega_i(t)+\epsilon_i}
with \eqn{Y_i(t)} the output at time \eqn{t} for the \eqn{i}th individual; \eqn{X_i(t)} the input predictors (fixed effects) at time \eqn{t} for the \eqn{i}th individual;
\eqn{Z_i(t)} are the random effects at time \eqn{t} for the \eqn{i}th individual; \eqn{\omega_i(t)} is a Brownian motion with volatility \eqn{\gamma^2=0.8} at time \eqn{t} for the \eqn{i}th individual; \eqn{\epsilon_i} is the residual error with
variance \eqn{\sigma^2=0.5}.
The data are simulated according to the simulations in low dimensional in the low dimensional scheme of the paper \url{doi:10.1177/0962280220946080}
}
\examples{
oldpar <- par()
oldopt <- options()
data <- DataLongGenerator(n=17, p=6,G=6) # Generate the data
# Let's see the output :
w <- which(data$id==1)
plot(data$time[w],data$Y[w],type="l",ylim=c(min(data$Y),max(data$Y)), col="grey")
for (i in unique(data$id)){
  w <- which(data$id==i)
  lines(data$time[w],data$Y[w], col='grey')
}
# Let's see the fixed effects predictors:
par(mfrow=c(2,3), mar=c(2,3,3,2))
for (i in 1:ncol(data$X)){
  w <- which(data$id==1)
  plot(data$time[w],data$X[w,i], col="grey",ylim=c(min(data$X[,i]),
  max(data$X[,i])),xlim=c(1,max(data$time)),main=latex2exp::TeX(paste0("$X^{(",i,")}$")))
  for (k in unique(data$id)){
    w <- which(data$id==k)
    lines(data$time[w],data$X[w,i], col="grey")
  }
}
par(oldpar)
options(oldopt)

}
