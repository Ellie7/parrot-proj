# Yellow-shouldered amazon life stage simulation analysis 
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

# making matrices using 'ysaFunc' function 
ysaFunc(yellow)
yellow
#make 10 matrices of random numbers.#for current code these aren't random and are drawn from ysa beta and normal distributed 
#vital rates
mat1 <- map(1:10, function(x) ysaFunc(yellow))

# use these matrices in mat 1 and get the eigen system for each....
mat2<-map(mat1, function(x) eigen(x)) 

#finding the first eigenvalue (finite rate of increase), lambda for each matrix 
mat3 <- map(mat2, function(x) {
  dom.pos <- which.max(x[["values"]])
  Re(x[["values"]][dom.pos])
})
#rs 
mat4 <- map(mat3, function(x) log(x)) 

##### alternatively 
analysis <- map(mat1, eigen.analysis) 
#rs
map(analysis, function(x) log(x$lambda1))

elasticities <- map(analysis, function(x) x$elasticities)

map(analysis, function(x) x$elasticities
  G1 <- x$elasticities[2,1] 
  G2 <- x$elasticities[3,2] 
  P2 <- x$elasticities[2,2] 
  P3 <- x$elasticities[3,3] 
  F3 <- x$elasticities[1,3]
  )
map(elasticities, function(x) 
  G1 <- x[2,1] 
  G2 <- x[3,2] 
  P2 <- x[2,2] 
  P3 <- x[3,3] 
  F3 <- x[1,3])

ggplot(elasticities, aes(x = vital rate, y = mean(elasticities))

       
n0 <- c(1, 1, 1)
stoch.projection(mat1, n0, tmax = 50, nreps = 10, prob = NULL,
                        nmax = NULL, sumweight = rep(1, length(n0)), verbose=FALSE)

    