#' (S)MERF algorithm
#'
#' (S)MERF is an adaptation of the random forest regression method to longitudinal data introduced by Hajjem et. al. (2014) <doi:10.1080/00949655.2012.741599>.
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
#' @param mtry [numeric]: Number of variables randomly sampled as candidates at each split. The default value is \code{p/3}.
#' @param ntree [numeric]: Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times. The default value is \code{ntree=500}.
#' @param time [vector]: Is the vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' @param sto [character]: Defines the covariance function of the stochastic process, can be either \code{"none"} for no stochastic process, \code{"BM"} for Brownian motion, \code{OrnUhl} for standard Ornstein-Uhlenbeck process, \code{BBridge} for Brownian Bridge, \code{fbm} for Fractional Brownian motion; can also be a function defined by the user.
#' @param delta [numeric]: The algorithm stops when the difference in log likelihood between two iterations is smaller than \code{delta}. The default value is set to O.O01
#'
#' @import randomForest
#' @import stats
#' @return A fitted (S)MERF model which is a list of the following elements: \itemize{
#' \item \code{forest:} Random forest obtained at the last iteration.
#' \item \code{random_effects :} Predictions of random effects for different trajectories.
#' \item \code{id_btilde:} Identifiers of individuals associated with the predictions \code{random_effects}.
#' \item \code{var_random_effects: } Estimation of the variance covariance matrix of random effects.
#' \item \code{sigma_sto: } Estimation of the volatility parameter of the stochastic process.
#' \item \code{sigma: } Estimation of the residual variance parameter.
#' \item \code{time: } The vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' \item \code{sto: } Stochastic process used in the model.
#' \item \code{Vraisemblance:} Log-likelihood of the different iterations.
#' \item \code{id: } Vector of the identifiers for the different trajectories.
#' \item \code{OOB: } OOB error of the fitted random forest at each iteration.
#' }
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' data <- DataLongGenerator(n=20) # Generate the data composed by n=20 individuals.
#' # Train a SMERF model on the generated data. Should take ~ 50 seconds
#' # The data are generated with a Brownian motion,
#' # so we use the parameter sto="BM" to specify a Brownian motion as stochastic process
#' smerf <- MERF(X=data$X,Y=data$Y,Z=data$Z,id=data$id,time=data$time,mtry=2,ntree=500,sto="BM")
#' smerf$forest # is the fitted random forest (obtained at the last iteration).
#' smerf$random_effects # are the predicted random effects for each individual.
#' smerf$omega # are the predicted stochastic processes.
#' plot(smerf$Vraisemblance) # evolution of the log-likelihood.
#' smerf$OOB # OOB error at each iteration.
#'
#'
MERF <- function(X,Y,id,Z,iter=100,mtry=ceiling(ncol(X)/3),ntree=500, time, sto, delta = 0.001){
  q <- dim(Z)[2]
  nind <- length(unique(id))
  btilde <- matrix(0,nind,q) #### Pour la ligne i, on a les effets al?atoires de l'individu i
  sigmahat <- 1 #### init
  Btilde <- diag(rep(1,q)) ### init
  epsilonhat <- rep(0,length(Y))
  id_btilde <- unique(id)
  Tiime <- sort(unique(time))
  omega <- rep(0,length(Y))
  sigma2 <- 1
  Vrai <- NULL
  inc <- 1
  OOB <- NULL

  if (class(sto)=="character"){
    if (sto=="fbm"){
      id_omega <- matrix(0,nind,length(unique(time)))
      for (i in 1:length(unique(id))){
        w <- which(id ==id_btilde[i])
        time11 <- time[w]
        where <- NULL
        for (j in 1:length(time11)){
          where <- c(where,which(Tiime==time11[j]))
        }
        id_omega[i,where] <- 1
      }
      omega <- matrix(0,nind,length(unique(time)))
      omega2 <- rep(0,length(Y))
      h <- opti.FBM(X,Y,id,Z,iter, mtry,ntree,time)
      for (i in 1:iter){
        ystar <- rep(0,length(Y))
        for (k in 1:nind){ #### on retrace les effets al?atoires
          indiv <- which(id==unique(id)[k])
          ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
        }


        forest <- randomForest(X,ystar,mtry=mtry,ntree=ntree, importance = TRUE) ### on construit l'arbre
        fhat <- predict(forest) #### pr?diction avec l'arbre
        OOB[i] <- forest$mse[ntree]
        for (k in 1:nind){ ### calcul des effets al?atoires par individu
          indiv <- which(id==unique(id)[k])
          K <- cov.fbm(time[indiv], h)
          V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))+ sigma2*K
          btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          omega[k,which(id_omega[k,]==1)] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          omega2[indiv] <- omega[k,which(id_omega[k,]==1)]
          epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[k,which(id_omega[k,]==1)]
        }
        sigm <- sigmahat
        B <- Btilde
        sigmahat <- sig.fbm(Y,sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,h) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
        Btilde  <- bay.fbm(btilde,Btilde,Z,id,sigm, time, sigma2,h) #### MAJ des param?tres de la variance des effets al?atoires.
        ### MAJ de la volatilit? du processus stochastique
        sigma2 <- gam_fbm(Y,sigm,id,Z,B,time,sigma2,omega,id_omega,h)
        Vrai <- c(Vrai, logV.fbm(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,h))
        if (i>1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
        if (inc < delta) {
          print(paste0("stopped after ", i, " iterations."))
          sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat,sigma_sto=sigma2, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time =time, Hurst=h, OOB =OOB, omega=omega2)
          class(sortie)<-"longituRF"
          return(sortie)
        }
      }
      sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id),sigma_sto=sigma2,omega=omega2, sigma_sto =sigma2, time = time, sto= sto, Hurst =h, id=id, Vraisemblance=Vrai, OOB =OOB)
      class(sortie) <- "longituRF"
      return(sortie)
    }


    if (sto=="exp"){
      id_omega <- matrix(0,nind,length(unique(time)))
      for (i in 1:length(unique(id))){
        w <- which(id ==id_btilde[i])
        time11 <- time[w]
        where <- NULL
        for (j in 1:length(time11)){
          where <- c(where,which(Tiime==time11[j]))
        }
        id_omega[i,where] <- 1
      }
      omega <- matrix(0,nind,length(unique(time)))
      omega2 <- rep(0,length(Y))
      alpha <- opti.exp(X,Y,id,Z,iter, mtry,ntree,time)
      for (i in 1:iter){
        ystar <- rep(0,length(Y))
        for (k in 1:nind){ #### on retrace les effets al?atoires
          indiv <- which(id==unique(id)[k])
          ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
        }

        forest <- randomForest(X,ystar,mtry=mtry,ntree=ntree, importance = TRUE) ### on construit l'arbre
        fhat <- predict(forest) #### pr?diction avec l'arbre
        OOB[i] <- forest$mse[ntree]
        for (k in 1:nind){ ### calcul des effets al?atoires par individu
          indiv <- which(id==unique(id)[k])
          K <- cov.exp(time[indiv], alpha)
          V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))+ sigma2*K
          btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          omega[k,which(id_omega[k,]==1)] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          omega2[indiv] <- omega[k,which(id_omega[k,]==1)]
          epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
        }
        sigm <- sigmahat
        B <- Btilde
        sigmahat <- sig.exp(Y,sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,alpha) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
        Btilde  <- bay.exp(btilde,Btilde,Z,id,sigm, time, sigma2,alpha) #### MAJ des param?tres de la variance des effets al?atoires.
        ### MAJ de la volatilit? du processus stochastique
        sigma2 <- gam_exp(Y,sigm,id,Z,B,time,sigma2,omega,id_omega,alpha)
        Vrai <- c(Vrai,logV.exp(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,alpha))
        if (i>1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
        if (inc < delta) {
          print(paste0("stopped after ", i, " iterations."))
          sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time=time, alpha = alpha, OOB =OOB, omega=omega2)
          class(sortie) <- "longituRF"
          return(sortie)
        }
      }
      sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), omega=omega2, sigma_sto =sigma2, time = time, sto= sto, alpha=alpha, id=id, Vraisemblance=Vrai, OOB =OOB)
      class(sortie) <- "longituRF"
      return(sortie)
    }

    if ( sto=="none"){
      for (i in 1:iter){
        ystar <- rep(NA,length(Y))
        for (k in 1:nind){ #### on retrace les effets al?atoires
          indiv <- which(id==unique(id)[k])
          ystar[indiv] <- Y[indiv]- Z[indiv,,drop=FALSE]%*%btilde[k,]
        }
        forest <- randomForest(X,ystar,mtry=mtry,ntree=ntree, importance = TRUE) ### on construit l'arbre
        fhat <- predict(forest)
        OOB[i] <- forest$mse[ntree]
        for (k in 1:nind){
          indiv <- which(id==unique(id)[k])
          V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))
          btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]
        }

        sigm <- sigmahat
        sigmahat <- sig(sigma = sigmahat,id = id, Z = Z, epsilon = epsilonhat,Btilde = Btilde)
        Btilde  <- bay(bhat = btilde,Bhat = Btilde,Z = Z,id = id,sigmahat = sigm)
        Vrai <- c(Vrai, logV(Y,fhat,Z,time,id,Btilde,0,sigmahat,sto))
        if (i>1) inc <- abs((Vrai[i-1]-Vrai[i])/Vrai[i-1])

        if (inc < delta) {
          print(paste0("stopped after ", i, " iterations."))
          sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time=time, OOB =OOB)
          class(sortie) <- "longituRF"
          return(sortie)
        }
      }
      sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance=Vrai,id=id, time=time, OOB =OOB)
      class(sortie) <- "longituRF"
      return(sortie)
    }
  }
  for (i in 1:iter){
    ystar <- rep(0,length(Y))
    for (k in 1:nind){
      indiv <- which(id==unique(id)[k])
      ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
    }
    forest <- randomForest(X,ystar,mtry=mtry,ntree=ntree, importance=TRUE)
    fhat <- predict(forest)
    OOB[i] <- forest$mse[ntree]
    for (k in 1:nind){
      indiv <- which(id==unique(id)[k])
      K <- sto_analysis(sto,time[indiv])
      V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))+ sigma2*K
      btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      omega[indiv] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
    }
    sigm <- sigmahat
    B <- Btilde
    sigmahat <- sig_sto(sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,sto)
    Btilde  <- bay_sto(btilde,Btilde,Z,id,sigm, time, sigma2,sto)
    sigma2 <- gam_sto(sigm,id,Z,B,time,sigma2,sto,omega)
    Vrai <- c(Vrai, logV(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,sto))
    if (i>1) inc <- abs((Vrai[i-1]-Vrai[i])/Vrai[i-1])
    if (inc < delta) {
      print(paste0("stopped after ", i, " iterations."))
      sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai,id=id, OOB =OOB)
      class(sortie) <- "longituRF"
      return(sortie)
    }
  }
  sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id),omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai,id=id, OOB =OOB)
  class(sortie) <- "longituRF"
  return(sortie)
}
