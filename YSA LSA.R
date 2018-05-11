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
map(mat1, eigen.analysis) 



