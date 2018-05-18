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

# making the function 
ysaFuncDD <- function (dataSource, n, threshold) 
{ 
  #ps
  p1a<- betaval((dataSource$pi[1]), (dataSource$piSD[1]), fx=runif(1)) 
  p1b<- betaval((dataSource$pi[2]), (dataSource$piSD[2]), fx=runif(1)) 
  p1c<- betaval((dataSource$pi[3]), (dataSource$piSD[3]), fx=runif(1))
  p2 <- betaval((dataSource$pi[4]), (dataSource$piSD[4]), fx=runif(1))
  p3 <- betaval((dataSource$pi[5]), (dataSource$piSD[5]), fx=runif(1))
  
  #f
  f3 <- rnorm(1, mean = (dataSource$f[5]), sd = (dataSource$fSD[5])) 
  
  # Pi <- ((1 - (pi^(di - 1)))/(1 - (pi^di)))*pi ------- equation for Pi's
  # Gi <- (pi^di*(1 - pi))/(1 - pi^di)           ------- equation for Gi's
  
  #d
  d1 <- dataSource$di[1] + dataSource$di[2] + dataSource$di[3]
  d2 <- dataSource$di[4]
  d3 <- dataSource$di[5]
  
  # this uses p1's defined above
  p1 <- (p1a*p1b*p1c) # this stage as the survival is from the multiplication of  p1a, p1b and p1c
  #add ps 
  
  # construct the matrix using defined parameters above
  matrix2 <- matrix(0, nrow = 3, ncol = 3)
  dimnames(matrix2) <- list(rownames(matrix2, do.NULL = FALSE, prefix = "row"),
                            colnames(matrix2, do.NULL = FALSE, prefix = "col"))
  matrix2[1,1] <- ((1 - (p1^(d1 - 1)))/(1 - (p1^d1)))*p1 
  matrix2[2,2] <- ((1 - (p2^(d2 - 1)))/(1 - (p2^d2)))*p2
  matrix2[3,3] <- ((1 - (p3^(d3 - 1)))/(1 - (p3^d3)))*p3
  
  #add f
  matrix2[1,3] <- f3
  
  #add gs 
  matrix2[2,1] <- (p1^d1*(1 - p1))/(1 - p1^d1) 
  matrix2[3,2] <- (p2^d2*(1 - p2))/(1 - p2^d2) 
  
  return(matrix2)
} 

# using the function 
threshold <- 500 # limiting number of nests 
n <- c(100, 100, 100, 100, 100) # population vector 

A <- ysaFuncDD(yellow)  

n <- A %*% n 


