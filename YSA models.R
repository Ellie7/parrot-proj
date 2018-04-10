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
library(reshape2)
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


#matrix from mean values
projmat <- ysameanFunc(yellow)

#project this matrix into the future 100 years 
#assumes starting pop of 2 individuals in each age class
N0 <- c(100,100,100)
nits <- 50 #how many years we want to project
tmp <- matrix(NA, nrow= nits, ncol=3) #number of columns = number of stages 
collect <- data.frame(tmp)
names(collect) <- c("E","J","A")
head(collect)

collect[1,] <- N0
head(collect)

#make a for loop 
for (a in 2:nits){
  collect[a,] <- projmat %*% t(collect[(a-1),]) #t stand for transpose (take the row put it on its side)
}
#tail gives the number of individuals in the first 6 years
head(collect)
#tail gives the number of individuals in the last 6 years
tail(collect)

#calculate population growth rate 
eigen(projmat)$values[1]
#0.9515817

#calculate the stable age distribution using eigen

Total <- mutate(collect, TotalPop = E+J+A, Time=1:nits)
Total
#graph with stage classes
long<-melt(Total,
           # this is the fixed variable
           id.vars = c("Time"),
           # to stack up
           measure.vars=c("E","J","A"),
           # label the stacked categories
           variable.name = "Stage",
           # label the stacked values
           value.name = "Abund")
# make graphs 
ggplot(long, aes(x =Time, y = log(Abund), group = Stage))+
  geom_line()+
  theme_classic() 
ggplot(long, aes(x =Time, y = Abund, group = Stage))+
  geom_line()+
  theme_classic()
