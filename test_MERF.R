devtools::load_all(".")
rm(list = ls())
library(dplyr)
load("train2018.Rdata")
train = train_2018 %>%
  select(
    'ID',
    'AGE',
    'SEXE',
    'Speed',
    'KM',
    'NO2',
    'O3',
    'PM10',
    'PM25',
    'T',
    'FF',
    'PSTAT',
    'U',
    'INS',
    'nbpartcl'
  )

moy <-apply(train[, 6:14], 2, mean)
sd <-apply(train[, 6:14], 2, sd)
train[, 6:14] <-as.data.frame(scale(train[, 6:14], center = moy, scale = sd))

#train[, 6:14] <-as.data.frame(scale(train[, 6:14]))
X_train = data.matrix(subset(train, select = -c(ID, Speed, AGE, SEXE, nbpartcl)))
Y_train = train$Speed
Z_train = as.matrix(cbind(rep(1, nrow(train)), 2 * runif(nrow(train))))
id_train = train$ID
time_train = train$KM

samples <- 1:100

id_trainb = id_train[samples]
X_trainb = X_train[samples,]
Y_trainb = Y_train[samples]
Z_trainb = Z_train[samples,]
time_trainb = time_train[samples]


X = X_trainb
Y = Y_trainb
Z = Z_trainb
id = id_trainb
time = time_trainb
sto = 'none'
iter = 2
mtry = 2
ntree = 500
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
  i <- 1
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
