# -*- coding: utf-8 -*-
rm(list=ls())
devtools::load_all(".")


load( "train2018.Rdata" )
head(train_2018)

train = train_2018[,c('ID','AGE','SEXE','Speed','KM','NO2','O3','PM10',
         'PM25','T','FF','PSTAT','U','INS','nbpartcl')]


# centrage et réduction

moy <- apply(train[,6:14], 2, mean) 
sd <- apply(train[,6:14], 2, sd)
train[,6: 14] <- as.data.frame(scale(train[,6: 14], center = moy, scale = sd))

train[,6:14] <- as.data.frame(scale(train[, 6: 14]))
X_train = data.matrix(subset(train, select = -c(ID,Speed,AGE,SEXE,nbpartcl)))
Y_train = train$Speed
Z_train = as.matrix(cbind(rep(1, nrow(train)), 2 * runif(nrow(train))))
id_train = train$ID
time_train = train$KM


# ci j'ai encore réduit la taille de l'échantllon pour tester le pgm

samples = 1:40
id=id_train[samples]
X=X_train[samples,]
Y=Y_train[samples]
Z=Z_train[samples,]
time=time_train[samples]


MERF(X = X,
     Y = Y,
     id = id,
     Z = Z,
     time = time,
     sto = 'none',
     iter = 2,
     mtry = 2,
     ntree = 100, 
     delta = 0.1)


