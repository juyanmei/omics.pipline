#-------------------------------------------------------------------------------#
# Copyright (c) 2017 Yanmei Ju (BGI-shenzhen). Allrights reserved.              #
# Created by Yanmei Ju (BGI-shenzhen) on 12/15/2017                             #
# This R program is using to class state2 in lasso model cvfold repeat times    #
# Args:                                                                         #
#   motu.prof: column is sample, row is motu                                    #
#   state.prof: column is state and row is sample                               #
#   cv.fold: number of cv.fold                                                  #
#   repeat.time: time of repeating cv.fold                                      #
#   prefix: the result prefix                                                   #
# output:                                                                       #
    # roc figure                                                                # 
    # rep.coef.nozero: final variable selected                                  #
# library(pROC)                                                                 #
# library(glmnet)                                                               #
#-------------------------------------------------------------------------------#

workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

# load package
library(glmnet)
source("../SZData/ROC.r")

#data preprocessing
motu.prof <- read.table("../SZData/mOTU.0.05.profile",header = 1,row.names = 1)
state.prof <- read.table("../SZData/state171.txt",header = 1,row.names = 1)
colnames(motu.prof) <- gsub("\\.","-",colnames(motu.prof))
state.prof <- state.prof[pmatch(colnames(motu.prof), rownames(state.prof)),]
motu.prof <- t(motu.prof)
cv.fold <- 10
repeat.time <- 10
prefix <- "lasso"

# function
cv.lasso <- function(x, y, cv.fold) {
  # cv lasso
  # Args:
  #  x: motu profile
  #  y: phenotype profile
  #  cv.fold: n fold CV
  
  # Returns:
  #  list1: lasso.pred; list2: lasso.coef
  N <- nrow(x)
  foldid <- sample(rep(seq(cv.fold), length = N))
  lasso.coef <- matrix(NA, ncol(x), cv.fold)
  lasso.pred <- matrix(NA, N, cv.fold)
  for (i in 1:cv.fold) {
    which = foldid == i
    if(is.matrix(y))y_sub = y[!which, ]else y_sub = y[!which]
    cv.out <- cv.glmnet(x[!which, ], y_sub, family = "binomial", type.measure="class", nfolds = 5,lambda = lams)
    lasso.pred[which, i] <- predict(cv.out$glmnet.fit, s = cv.out$lambda.min, newx = x[which, ], type="response")
    lasso.coef[, i] <- predict(cv.out$glmnet.fit, type='coefficients', s=cv.out$lambda.min)[1:ncol(x),]
  }
  rownames(lasso.pred) <- rownames(x)
  rownames(lasso.coef) <- colnames(x)
  list(predicted = lasso.pred, coef = lasso.coef)
}

#classification
lams =10^ seq (10, -2, length =500)
set.seed(1)
y <- state.prof[,1]
a <- replicate(repeat.time, cv.lasso(motu.prof, y, cv.fold))

# modify result
rep.pred <- NULL
rep.coef <- NULL
for (i in 1:repeat.time) {
  rep.pred <- cbind(rep.pred, a[, i]$predicted)
  rep.coef <- cbind(rep.coef, a[, i]$coef)
}
name <- matrix(NA, cv.fold, repeat.time)
for (i in 1:repeat.time) {
  for (j in 1:cv.fold) {
    name[j, i] <- paste("cv_fold", j, "rep", i, sep = "_")
  }
}
name <- matrix(name, 1, cv.time * repeat.time)
colnames(rep.coef) <- name

# nonzero coefficient number in samples at least 5: nozero.num
# nonzero coefficient in at least 50% of the LASSO models: nozero.var
nozero.var <- apply(rep.coef, 1, function(x) length(which(x != 0 )))
nozero.num <- apply(rep.coef, 2, function(x) length(which(x != 0 )))
rep.coef.nozero <- rep.coef[which(nozero.var >= (0.5 * cv.fold * repeat.time)), which(nozero.num >= 5)]
coef.name <- paste(prefix, "txt", sep = ".")
write.table(rep.coef.nozero, coef.name, quote = F, row.names = T, sep = "\t")

# ROC
final.pred <- apply(rep.pred, 1, function(x) mean(x, na.rm = T))
pdf.name <- paste(prefix, "pdf", sep = ".")
pdf(pdf.name)
plot_roc(y, final.pred)
dev.off()


