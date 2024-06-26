% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/REEMtree.R
\name{REEMtree}
\alias{REEMtree}
\title{(S)REEMtree algorithm}
\usage{
REEMtree(X, Y, id, Z, iter = 10, time, sto, delta = 0.001)
}
\arguments{
\item{X}{\link{matrix}: A \code{N}x\code{p} matrix containing the \code{p} predictors of the fixed effects, column codes for a predictor.}

\item{Y}{\link{vector}: A vector containing the output trajectories.}

\item{id}{\link{vector}: Is the vector of the identifiers for the different trajectories.}

\item{Z}{\link{matrix}: A \code{N}x\code{q} matrix containing the \code{q} predictor of the random effects.}

\item{iter}{\link{numeric}: Maximal number of iterations of the algorithm. The default is set to \code{iter=100}}

\item{time}{\link{vector}: Is the vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.}

\item{sto}{\link{character}: Defines the covariance function of the stochastic process, can be either \code{"none"} for no stochastic process, \code{"BM"} for Brownian motion, \code{OrnUhl} for standard Ornstein-Uhlenbeck process, \code{BBridge} for Brownian Bridge, \code{fbm} for Fractional Brownian motion; can also be a function defined by the user.}

\item{delta}{\link{numeric}: The algorithm stops when the difference in log likelihood between two iterations is smaller than \code{delta}. The default value is set to O.O01}
}
\value{
A fitted (S)MERF model which is a list of the following elements: \itemize{
\item \code{forest:} Tree obtained at the last iteration.
\item \code{random_effects :} Predictions of random effects for different trajectories.
\item \code{id_btilde:} Identifiers of individuals associated with the predictions \code{random_effects}.
\item \code{var_random_effects: } Estimation of the variance covariance matrix of random effects.
\item \code{sigma_sto: } Estimation of the volatility parameter of the stochastic process.
\item \code{sigma: } Estimation of the residual variance parameter.
\item \code{time: } The vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
\item \code{sto: } Stochastic process used in the model.
\item \code{Vraisemblance:} Log-likelihood of the different iterations.
\item \code{id: } Vector of the identifiers for the different trajectories.
}
}
\description{
(S)REEMtree is an adaptation of the random forest regression method to longitudinal data introduced by Sela and Simonoff. (2012) \url{doi:10.1007/s10994-011-5258-3}.
The algorithm will estimate the parameters of the following semi-parametric stochastic mixed-effects model: \deqn{Y_i(t)=f(X_i(t))+Z_i(t)\beta_i + \omega_i(t)+\epsilon_i}
with \eqn{Y_i(t)} the output at time \eqn{t} for the \eqn{i}th individual; \eqn{X_i(t)} the input predictors (fixed effects) at time \eqn{t} for the \eqn{i}th individual;
\eqn{Z_i(t)} are the random effects at time \eqn{t} for the \eqn{i}th individual; \eqn{\omega_i(t)} is the stochastic process at time \eqn{t} for the \eqn{i}th individual
which model the serial correlations of the output measurements; \eqn{\epsilon_i} is the residual error.
}
\examples{
set.seed(123)
data <- DataLongGenerator(n=20) # Generate the data composed by n=20 individuals.
# Train a SREEMtree model on the generated data.
# The data are generated with a Brownian motion,
# so we use the parameter sto="BM" to specify a Brownian motion as stochastic process
X.fixed.effects <- as.data.frame(data$X)
sreemt <- REEMtree(X=X.fixed.effects,Y=data$Y,Z=data$Z,id=data$id,time=data$time,
sto="BM", delta=0.0001)
sreemt$forest # is the fitted random forest (obtained at the last iteration).
sreemt$random_effects # are the predicted random effects for each individual.
sreemt$omega # are the predicted stochastic processes.
plot(sreemt$Vraisemblance) #evolution of the log-likelihood.

}
