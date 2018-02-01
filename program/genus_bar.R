#-------------------------------------------------------------------------------#
# Copyright (c) 2018 Yanmei Ju (BGI-shenzhen). Allrights reserved.              #
# Created by Yanmei Ju (BGI-shenzhen) on 02/01/2018                             #
# This R program is using to draw taxo relative abundance bar plot              #
# Attention: top.num should less than 12, this progrom need to updata           #
# Args:                                                                         #
#   taxo.prof: row is sample, column is taxo and the first column is no rank    #
#   state.prof: column is state and row is sample                               #
#   time.point: the time point(int: 1-5)                                        #
#   prefix: output name                                                         #
#   top.num: the top top.num(int) taxo number                                   #
#   norank: T/F(T for figure includes no rank, F for figure not include no rank)#
# output:                                                                       #
#   out: figure of top.num relative abundance                                   #   
# library(ggplot2)                                                              #
# library(reshape2)
#-------------------------------------------------------------------------------#

# load dir
workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

# load package
library(ggplot2)
library(reshape2)


# load data
prefix <- "genus.16S.reabun"
time.point <- 5
top.num <- 10
norank <- F
taxo.prof <- read.table("genus.16s.abun.txt", header = 1, row.names = 1)
state.prof <- read.table("group.time.txt", header = 1, row.names = 1)

# function
groupSplit <- function(x, y, time, group) {
  # x: state, row is sample
  # y: taxo.prof, col is taxo; row is sample
  state.x <- x[which(x$time == time & x$group == group), ]
  y[pmatch(rownames(state.x), rownames(y)), ]
}
top <- function(x, top) {
  x.sort <- sort(x, decreasing = T)
  x.sort[1:top]
}

################   data preprocessing    ###################
# choose time.point profile
taxo.prof.case <- groupSplit(state.prof, taxo.prof, time.point, 'case')
taxo.prof.control <- groupSplit(state.prof, taxo.prof, time.point, 'control')
taxo.prof.c <- rbind(taxo.prof.case, taxo.prof.control)

# choose top top.num taxo
name <- NULL
for (i in 1:nrow(taxo.prof.c)) {
  x <- taxo.prof.c[i, ]
  x <- x[which(x != 0)]
  x.top <- top(x, top.num)
  name <- c(name, names(x.top))
}
name.f <- top(table(name), top.num + 1)

# modify result for drawing figure
taxo.final <- taxo.prof.c[, pmatch(names(name.f), colnames(taxo.prof.c))]
if(norank) {
  other <- apply(taxo.final, 1, function(x) 1-sum(x))
  total <- cbind(rownames(taxo.final), taxo.final, other)
} else if(norank == F) {
  taxo.final <- taxo.final[, -which(colnames(taxo.final) == 'No_Rank')]
  other <- apply(taxo.final, 1, function(x) 1-sum(x))
  total <- cbind(rownames(taxo.final),taxo.final, other)
}
total <- melt(total)
colnames(total) <- c("id", "bac", "val")

# draw the figure
fig.name <- paste(prefix, time.point, "pdf", sep = ".")
ggplot(total, aes(x = id, y = val, fill = bac)) +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_brewer(palette = "Set3") +
  ylab("relative abundance") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
ggsave(fig.name, width = 15, height = 12, units = "cm")











