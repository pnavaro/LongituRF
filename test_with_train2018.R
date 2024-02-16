# -*- coding: utf-8 -*-
rm(list=ls())
# Import and install library
library(dplyr)
# install.packages("LongituRF")
library(LongituRF)


load( "train2018.Rdata" )
head(train_2018)

train = train_2018 %>%
  select('ID','AGE','SEXE','Speed','KM','NO2','O3','PM10',
         'PM25','T','FF','PSTAT','U','INS','nbpartcl')


# centrage et réduction

moy <- apply(train[,6:14], 2, mean) #il y avait une erreur ici avec l'oubli de la virgule
sd <- apply(train[,6:14], 2, sd) #il y avait une erreur ici avec l'oubli de la virgule
train[,6: 14] <- as.data.frame(scale(train[,6: 14], center = moy, scale = sd))

train[,6:14] <- as.data.frame(scale(train[, 6: 14]))#il y avait une erreur ici avec l'oubli de la virgule
X_train = data.matrix(subset(train, select = -c(ID,Speed,AGE,SEXE,nbpartcl)))
Y_train = train$Speed
Z_train = as.matrix(cbind(rep(1, nrow(train)), 2 * runif(nrow(train))))
id_train = train$ID
time_train = train$KM


# ci j'ai encore réduit la taille de l'échantllon pour tester le pgm

id_trainb=id_train[1:40]
X_trainb=X_train[1:40,]
Y_trainb=Y_train[1:40]
Z_trainb=Z_train[1:40,]
time_trainb=time_train[1:40]


MERF(X = X_trainb,
     Y = Y_trainb,
     id = id_trainb,
     Z = Z_trainb,
     time = time_trainb,
     sto = 'none',
     iter=2,
     mtry = 2,
     ntree = 100,delta=0.1)
