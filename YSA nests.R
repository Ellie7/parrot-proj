####### YSA nest limitation modelling 
# apply a threshold based density dependence to the fecundity term 
# Once the number of breeding pairs reaches the number of nests available, the total reproductive output can’t increase any more.
# The way we would do that is to replace the fixed fecundity term by a function that models the per-capita recruitment based on that
# threshold. That means that you have to re-calculate the matrix at each time step based on the current population vector, so we can’t
# use the same type of analysis as we’ve used for the fixed matrix model.
# So, we can’t calculate the max eigen value for the matrix because the matrix varies with population density. But we can look at long 
# term behaviour by running simulation at looking at what we get - we could look at population structure and population size (which 
# will be limited).
# Then we can look at how different thresholds influence these things and what that suggests will happen if the number of nests is 
# being reduced through poaching.

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
source(file = "YSA life history data.R")

# N > M  -> (M*F)/N
# N <= M ->  F

threshold <- 500 # limiting number of nests 
n <- c(100, 100, 100, 100, 100) # population vector 
  
ysaFunc <- function(dataSource)
ysaFuncDD <- function (dataSource, n, threshold)
  
A <- ysaFuncDD(yellow)  

n <- A %*% n


