devtools::load_all(".")

data <- DataLongGenerator(n=20) # Generate the data composed by n=20 individuals.

X <- data$X
Y <- data$Y
Z <- data$Z
id <- data$id
time <- data$time
mtry <- 2
ntree <- 500

sto = 'none'
iter = 2
mtry = 2
delta = 0.1

q <- dim(Z)[2]
nind <- length(unique(id))
btilde <-matrix(0, nind, q)
sigmahat <- 1
Btilde <- diag(rep(1, q))
epsilonhat <- rep(0, length(Y))
id_btilde <- unique(id)
Tiime <- sort(unique(time))
omega <- rep(0, length(Y))
sigma2 <- 1
Vrai <- NULL
inc <- 1
OOB <- NULL

for (i in 1:iter) {
  ystar <- rep(NA, length(Y))
  print(paste(" i = ", i))
  for (k in 1:nind) {
    indiv <- which(id == unique(id)[k])
    ystar[indiv] <-
      Y[indiv] - Z[indiv, , drop = FALSE] %*% btilde[k,]
  }
  print("ystar=")
  print(ystar)
  forest <-
    randomForest(X,
                 ystar,
                 mtry = mtry,
                 ntree = ntree,
                 importance = TRUE)
  fhat <- predict(forest)
  print("fhat=")
  print(fhat)
  OOB[i] <- forest$mse[ntree]
  for (k in 1:nind) {
    indiv <- which(id == unique(id)[k])
    V <-
      Z[indiv, , drop = FALSE] %*% Btilde %*% t(Z[indiv, , drop = FALSE]) + diag(as.numeric(sigmahat), length(indiv), length(indiv))
    btilde[k,] <-
      Btilde %*% t(Z[indiv, , drop = FALSE]) %*% solve(V) %*% (Y[indiv] - fhat[indiv])
    epsilonhat[indiv] <-
      Y[indiv] - fhat[indiv] - Z[indiv, , drop = FALSE] %*% btilde[k,]
  }

  sigm <- sigmahat
  sigmahat <-
    sig(
      sigma = sigmahat,
      id = id,
      Z = Z,
      epsilon = epsilonhat,
      Btilde = Btilde
    )
  Btilde  <-
    bay(
      bhat = btilde,
      Bhat = Btilde,
      Z = Z,
      id = id,
      sigmahat = sigm
    )
  Vrai <-
    c(Vrai, logV(Y, fhat, Z, time, id, Btilde, 0, sigmahat, sto))
  print(paste0("Vrai = ", Vrai[i]))

  if (i > 1) {

    inc <- abs((Vrai[i - 1] - Vrai[i]) / Vrai[i - 1])

    if (inc < delta) {
      print(paste0("stopped after ", i, " iterations."))
      sortie <-
        list(
          forest = forest,
          random_effects = btilde,
          var_random_effects = Btilde,
          sigma = sigmahat,
          id_btilde = unique(id),
          sto = sto,
          vraisemblance = Vrai,
          id = id,
          time = time,
          OOB = OOB
        )
      class(sortie) <- "longituRF"
      print(sortie)
      break
    }
  }
}
