#-------------------------------------------------------------------------------#
# Copyright (c) 2018 Yanmei Ju (BGI-shenzhen). Allrights reserved.              #
# Created by Yanmei Ju (BGI-shenzhen) on 01/09/2018                             #
# This R program is using to analyse phenotype significant                      #
# Args:                                                                         #
#   motu.prof: column is sample, row is motu                                    #
#   phe.prof: column is phe and row is sample                                       
#   ps: permutation time
#   prefix: the result prefix                                                   #
#   phe.type: the type of phenotype, such as numeric, factor                    #
# output:                                                                       #
#   out: list such as R2, pvalue                                                #   
# library(vegan)                                                                 #
#-------------------------------------------------------------------------------#

# load dir
workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

# load packages
library(vegan)

# load data
ps <- 999
prefix <- "envfit.txt"
motu.prof <- read.table("../SZData/mOTU.0.05.profile", header = 1, row.names = 1)
phe.type <- read.table("phe_type_1213.txt", stringsAsFactors = F)
phe.prof <- read.table("total_phe_1213.txt", header = 1, row.names = 1, colClasses = c(NA, phe.type[, 1]))
colnames(motu.prof) <- gsub("\\.", "-", colnames(motu.prof))
phe.prof <- phe.prof[pmatch(colnames(motu.prof), rownames(phe.prof)), ]

# data preprocessing
motu.prof <- t(motu.prof)
phe.cn <- colnames(phe.prof)

# Ordination
ord <- metaMDS(motu.prof)
#ord <- cca(motu.prof)

out.cn <- list("Phenotype", "SampleNum", "type", "R2", "Pvalue")
write.table(out.cn, prefix, quote = F, sep = "\t", col.names = F, row.names = F, append = F)
for (i in 1:ncol(phe.prof)) {
  phe <- phe.prof[, i]
  flag <- which(!is.na(phe))
  len <- length(flag)
  set.seed(0)
  res <- envfit(ord ~ phe, permutations = ps, na.rm = T)
  if (is.factor(phe)) {
    type <- "factor"
    r2 <- res$factors$r
    pvalue <- res$factors$pvals
  } else {
    type <- "numeric"
    r2 <- res$vectors$r
    pvalue <- res$vectors$pvals
  }
  out <- list(phe.cn[i], len, type, r2, pvalue)
  write.table(out, prefix, quote = F, sep = "\t", col.names = F, row.names = F, append = T)
}





