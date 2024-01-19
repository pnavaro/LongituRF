#' Title
#'
#' @param X
#' @param Y
#' @param id
#' @param Z
#' @param iter
#' @param mtry
#' @param ntree
#' @param time
#'
#' @import stats
#'
#' @keywords internal
opti.FBMreem <- function(X,Y,id,Z,iter,mtry,ntree,time){
  print("Do you want to enter a set for the Hurst parameter ? (1/0)")
  resp <- scan(nmax=1)
  if (resp ==1){
    print("please enter your ensemble (vector):")
    H <- scan()
    opti <- NULL}

  if (resp==0) {H <- seq(0.1,0.9,0.1)}
  opti <- NULL
  for (h in H){
    q <- dim(Z)[2]
    nind <- length(unique(id))
    btilde <- matrix(0,nind,q) #### Pour la ligne i, on a les effets al?atoires de l'individu i
    sigmahat <- 1 #### init
    Btilde <- diag(rep(1,q)) ### init
    epsilonhat <- 0
    id_btilde <- unique(id)
    Tiime <- sort(unique(time))
    omega <- matrix(0,nind,length(unique(time)))
    mu = rep(0,length(id_btilde))
    sigma2 <- 1
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
    for (i in 1:iter){
      ystar <- rep(0,length(Y))
      for (k in 1:nind){ #### on retrace les effets al?atoires
        indiv <- which(id==unique(id)[k])
        ystar[indiv] <- Y[indiv]- Z[indiv,]%*%btilde[k,]- omega[k,which(id_omega[k,]==1)]
      }
      forest <- randomForest(as.data.frame(X), ystar,mtry=mtry,ntree=ntree)
      f1 <- predict(forest, as.data.frame(X), nodes=TRUE)
      trees <- attributes(f1)
      K <- ntree

      for (k in 1:K){
        Phi <- matrix(0,length(Y),length(unique(trees$nodes[,k])))
        indii <- which(forest$forest$nodestatus[,k]==-1)
        for (l in 1:dim(Phi)[2]){
          w <- which(trees$nodes[,k]==indii[l])
          Phi[w,l] <- 1
        }
        beta <- Moy_fbm(id,Btilde,sigmahat,Phi,Y,Z, h ,time, sigma2)
        forest$forest$nodepred[indii,k] <- beta
      }

      fhat <- predict(forest,as.data.frame(X))

      for (k in 1:nind){ ### calcul des effets al?atoires par individu
        indiv <- which(id==unique(id)[k])
        K <- cov.fbm(time[indiv],h)
        V <- Z[indiv,]%*%Btilde%*%t(Z[indiv,])+diag(rep(sigmahat,length(Y[indiv])))+ sigma2*K
        btilde[k,] <- Btilde%*%t(Z[indiv,])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      }
      #### pr?diction du processus stochastique:
      for (k in 1:length(id_btilde)){
        indiv <- which(id==unique(id)[k])
        K <- cov.fbm(time[indiv], h)
        V <- Z[indiv,]%*%Btilde%*%t(Z[indiv,])+diag(rep(sigmahat,length(Y[indiv])))+sigma2*K
        omega[k,which(id_omega[k,]==1)] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      }

      for (k in 1:nind){
        indiv <- which(id==unique(id)[k])
        epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,]%*%btilde[k,]- omega[k,which(id_omega[k,]==1)]
      }
      sigm <- sigmahat
      B <- Btilde
      sigmahat <- sig.fbm(Y,sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,h) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
      Btilde  <- bay.fbm(btilde,Btilde,Z,id,sigm, time, sigma2,h) #### MAJ des param?tres de la variance des effets al?atoires.
      ### MAJ de la volatilit? du processus stochastique
      sigma2 <- gam_fbm(Y,sigm,id,Z,B,time,sigma2,omega,id_omega, h)
    }
    opti <- c(opti,sigmahat)
  }
  return(H[which.min(opti)])
}
