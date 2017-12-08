#-------------------------------------------------------------------------------------#
# Copyright (c) 2017 Yanmei Ju (BGI-shenzhen). Allrights reserved.                    #
# Created by Yanmei Ju (BGI-shenzhen) on 12/8/2017                                    #
# This R program is using to calculate and draw db-RDA                                #
# Args:                                                                               #
#   genus: genus profile, column is sample and row is genus                           #
#   sample.phe: phenotypes profile, column is phenotype and row is sample             #
#¡¡ sample.state: the state of sample, column is state(factor) and row is sample      #
#   prefix: output file prefix                                                        #
#   top: show the largest CAP genus                                                   #
# output: figure of db-RDA                                                            # 
# require(vegan)                                                                      #
# require(ggplot2)                                                                    #
#-------------------------------------------------------------------------------------#

# set path
workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

# load library
library(vegan)
library(ggplot2)

# load data
genus <- t(read.table("genus.profile", header = T, row.names = 1))
sample.phe <- read.table("mccb.txt", header = T, row.names = 1)
sample.state <- read.table("../state171.txt", header = T, row.names = 1)
prefix <- "dbRDA_1208"
top <- 5

# match phe and genus profile
rownames(sample.phe) <- sub("\\-",".",rownames(sample.phe))
genus.phe <- genus[pmatch(rownames(sample.phe), rownames(genus)),]
sample.phe <- as.matrix(sample.phe)

# remove na
id <- rowSums(is.na(sample.phe)) == 0
sample.phe <- sample.phe[id, ]
genus.phe <- genus.phe[id, ]
sample.state <- sample.state[id, ]

# method2 dbrda
dbrda2 <- capscale(genus.phe ~ sample.phe, distance = "bray")
anova(dbrda2, step = 1000, perm.max = 1000)
dbrda2.spe <- scores(dbrda2, display = "sp")
dbrda2.samp <- scores(dbrda2, display = "sites")
dbrda2.env <- scores(dbrda2, display = "bp")

dbrda2.spe.order <- dbrda2.spe[order(abs(dbrda2.spe[,1]), abs(dbrda2.spe[,2]), decreasing = T),]
dbrda2.spe.top <- dbrda2.spe.order[1:top, ] 

# figure
fig.name <- paste(prefix, "pdf", sep = ".")
total <- cbind(dbrda2.samp, sample.state)
colnames(total) <- c("CAP1", "CAP2", "state2", "state3")
total$state2 <- as.factor(total$state2)
ggplot(total, aes(CAP1, CAP2, colour = state2)) +
  geom_point() + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  annotate("text", x = dbrda2.spe.top[,1], y = dbrda2.spe.top[,2], 
               label = "*") +
  annotate("text", x = dbrda2.spe.top[,1]+1, y = dbrda2.spe.top[,2], 
           label = rownames(dbrda2.spe.top))
ggsave(fig.name, width = 6, height = 5)



