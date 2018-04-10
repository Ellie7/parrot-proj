# Yellow-shouldered amazon demographic modelling  
rm(list=ls())
library(dplyr)
library(ggplot2) 
library(knitr)
library(primer)
library(Rmisc)
library(agricolae)
library(popbio)
library(MASS) 
library(tidyverse) 
source(file = "ysa functions.R")
source(file = "ysa demography.R")

############################################################## EFFECTS OF PER CAPITA GROWTH RATE (effects of different lambda)
mat2
#finding the first eigenvalue (finite rate of increase)
B <- ysaFunc(yellow) #different matrices generated each time for now 
C <- ysaFunc(yellow)
D <- ysaFunc(yellow)
E <- ysaFunc(yellow)
F <- ysaFunc(yellow)

eigs.B <- eigen(B)
eigs.C <- eigen(C)
eigs.D <- eigen(D)
eigs.E <- eigen(E)
eigs.F <- eigen(F)

bdom.pos <- which.max(eigs.B[["values"]])
cdom.pos <- which.max(eigs.C[["values"]])
ddom.pos <- which.max(eigs.D[["values"]])
edom.pos <- which.max(eigs.E[["values"]])
fdom.pos <- which.max(eigs.F[["values"]])

Lb <- Re(eigs.B[["values"]][bdom.pos])
Lc <- Re(eigs.C[["values"]][cdom.pos])
Ld <- Re(eigs.D[["values"]][ddom.pos])
Le <- Re(eigs.E[["values"]][edom.pos])
Lf <- Re(eigs.F[["values"]][fdom.pos])

N0 <- 100
time <- 0:10
lambdas <- c(Lb, Lc, Ld, Le, Lf, LFi, L1ai, L1bi, L1ci, L2i, L3i)
#use sapply to apply geometrix growth function to each lambda, x stands for each lambda, which the funcion uses to 
# calculate population size
N.all <-sapply(lambdas, function(x) N0*x^time)
matplot(time, N.all, xlab="Years", ylab="N", pch = 1:3)