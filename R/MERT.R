#' (S)MERT algorithm
#'
#'
#' (S)MERT is an adaptation of the random forest regression method to longitudinal data introduced by Hajjem et. al. (2011) <doi:10.1016/j.spl.2010.12.003>.
#' The model has been improved by Capitaine et. al. (2020) <doi:10.1177/0962280220946080> with the addition of a stochastic process.
#' The algorithm will estimate the parameters of the following semi-parametric stochastic mixed-effects model: \deqn{Y_i(t)=f(X_i(t))+Z_i(t)\beta_i + \omega_i(t)+\epsilon_i}
#' with \eqn{Y_i(t)} the output at time \eqn{t} for the \eqn{i}th individual; \eqn{X_i(t)} the input predictors (fixed effects) at time \eqn{t} for the \eqn{i}th individual;
#' \eqn{Z_i(t)} are the random effects at time \eqn{t} for the \eqn{i}th individual; \eqn{\omega_i(t)} is the stochastic process at time \eqn{t} for the \eqn{i}th individual
#'  which model the serial correlations of the output measurements; \eqn{\epsilon_i} is the residual error.
#'
#' @param X [matrix]: A \code{N}x\code{p} matrix containing the \code{p} predictors of the fixed effects, column codes for a predictor.
#' @param Y [vector]: A vector containing the output trajectories.
#' @param id [vector]: Is the vector of the identifiers for the different trajectories.
#' @param Z [matrix]: A \code{N}x\code{q} matrix containing the \code{q} predictor of the random effects.
#' @param iter [numeric]: Maximal number of iterations of the algorithm. The default is set to \code{iter=100}
#' @param time [vector]: Is the vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' @param sto [character]: Defines the covariance function of the stochastic process, can be either \code{"none"} for no stochastic process, \code{"BM"} for Brownian motion, \code{OrnUhl} for standard Ornstein-Uhlenbeck process, \code{BBridge} for Brownian Bridge, \code{fbm} for Fractional Brownian motion; can also be a function defined by the user.
#' @param delta [numeric]: The algorithm stops when the difference in log likelihood between two iterations is smaller than \code{delta}. The default value is set to O.O01
#'
#'
#' @import stats
#' @import rpart
#'
#'
#'
#' @return A fitted (S)MERF model which is a list of the following elements: \itemize{
#' \item \code{forest:} Tree obtained at the last iteration.
#' \item \code{random_effects :} Predictions of random effects for different trajectories.
#' \item \code{id_btilde:} Identifiers of individuals associated with the predictions \code{random_effects}.
#' \item \code{var_random_effects: } Estimation of the variance covariance matrix of random effects.
#' \item \code{sigma_sto: } Estimation of the volatility parameter of the stochastic process.
#' \item \code{sigma: } Estimation of the residual variance parameter.
#' \item \code{time: } The vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' \item \code{sto: } Stochastic process used in the model.
#' \item \code{Vraisemblance:} Log-likelihood of the different iterations.
#' \item \code{id: } Vector of the identifiers for the different trajectories.
#' }
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' data <- DataLongGenerator(n=20) # Generate the data composed by n=20 individuals.
#' # Train a SMERF model on the generated data. Should take ~ 50 secondes
#' # The data are generated with a Brownian motion,
#' # so we use the parameter sto="BM" to specify a Brownian motion as stochastic process
#' smert <- MERF(X=data$X,Y=data$Y,Z=data$Z,id=data$id,time=data$time,sto="BM")
#' smert$forest # is the fitted random forest (obtained at the last iteration).
#' smert$random_effects # are the predicted random effects for each individual.
#' smert$omega # are the predicted stochastic processes.
#' plot(smert$Vraisemblance) #evolution of the log-likelihood.
#'
#'
MERT <- function(X,Y,id,Z,iter=100,time, sto, delta = 0.001){
  q <- dim(Z)[2]
  nind <- length(unique(id))
  btilde <- matrix(0,nind,q) #### Pour la ligne i, on a les effets al?atoires de l'individu i
  sigmahat <- 1 #### init
  Btilde <- diag(rep(1,q)) ### init
  epsilonhat <- 0
  id_btilde <- unique(id)
  Tiime <- sort(unique(time))
  omega <- rep(0,length(Y))
  sigma2 <- 1
  inc <- 1
  Vrai <- NULL
  id_omega=sto

  if (class(sto)=="character"){
    if ( sto=="none"){
      for (i in 1:iter){
        ystar <- rep(0,length(Y))
        for (k in 1:nind){ #### on retrace les effets al?atoires
          indiv <- which(id==unique(id)[k])
          ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]
        }

        tree <- rpart(ystar~.,as.data.frame(X)) ### on construit l'arbre
        fhat <- predict(tree, as.data.frame(X)) #### pr?diction avec l'arbre
        for (k in 1:nind){ ### calcul des effets al?atoires par individu
          indiv <- which(id==unique(id)[k])
          V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))
          btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]
        }
        sigm <- sigmahat
        sigmahat <- sig(sigmahat,id, Z, epsilonhat, Btilde) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
        Btilde  <- bay(btilde,Btilde,Z,id,sigm) #### MAJ des param?tres de la variance des effets al?atoires.
        Vrai <- c(Vrai, logV(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,0,sigmahat,"none"))
        if (i>1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
        if (inc< delta) {
          print(paste0("stopped after ", i, " iterations."))
          sortie <- list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, Vraisemblance = Vrai, id =id, time=time)
          class(sortie) <- "longituRF"
          return(sortie)
        }
      }
      sortie <- list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, Vraisemblance=Vrai, id=id, time=time)
      class(sortie) <- "longituRF"
      return(sortie)
    }
  }
  for (i in 1:iter){
    ystar <- rep(0,length(Y))
    for (k in 1:nind){ #### on retrace les effets al?atoires
      indiv <- which(id==unique(id)[k])
      ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
    }

    tree <- rpart(ystar~.,X) ### on construit l'arbre
    fhat <- predict(tree, X) #### pr?diction avec l'arbre
    for (k in 1:nind){ ### calcul des effets al?atoires par individu
      indiv <- which(id==unique(id)[k])
      K <- sto_analysis(sto,time[indiv])
      V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))+ sigma2*K
      btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      omega[indiv] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv]-Z[indiv,,drop=FALSE]%*%btilde[k,])
      epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
    }
    #### pr?diction du processus stochastique:
    sigm <- sigmahat
    B <- Btilde
    sigmahat <- sig_sto(sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,sto) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
    Btilde  <- bay_sto(btilde,Btilde,Z,id,sigm, time, sigma2,sto) #### MAJ des param?tres de la variance des effets al?atoires.
    ### MAJ de la volatilit? du processus stochastique
    sigma2 <- gam_sto(sigm,id,Z,B,time,sigma2,sto,omega)
    Vrai <- c(Vrai, logV(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,sto))
    if (i>1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
    if (inc< delta) {
      print(paste0("stopped after ", i, " iterations."))
      sortie <- list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat,id_omega=id_omega, id_btilde=unique(id),omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai, id = id)
      class(sortie) <- "longituRF"
      return(sortie)
    }
  }
  sortie <- list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id),id_omega=id_omega,omega=omega, sigma_sto =sigma2, time = time, sto= sto, Vraisemblance=Vrai, id=id)
  class(sortie) <- "longituRF"
  return(sortie)
}
