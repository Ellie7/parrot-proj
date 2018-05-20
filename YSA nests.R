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
# N > M  -> (M*F)/N
# N <= M ->  F

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
source(file = "YSA demography.R")
source(file = "YSA life history data.R")

# using ysaFuncDD function 

ysaFuncDD

# using the function 
threshold <- 500 # limiting number of nests 
n <- c(100, 100, 100) # population vector 

A <- ysaFuncDD(yellow, n, threshold, stochastic = FALSE)  

n <- A %*% n 

## Simulation with density dependence
# used to find out the carrying capacity for a given threshold number of nests.
sim_len <- 200
N <- vector(mode = "list", length = sim_len)
thres = 150
n <- c(50, 50, 50)
for (i in 1:sim_len) {
  A <- ysaFuncDD(yellow, n, thres)
  n <- A %*% n
  N[[i]] <- n
}
map(N, sum) %>% unlist %>% plot


#use an implementation of the function described above to simulate population trajectories under a range of poaching 
# probabilities: .01%, 0.1%, 0.5%, 1%
# used to simulate the population size under different levels of poaching related nest site loss at say 50 or 100 years in the 
#future.
# 0.01% poaching 
## Simulation with losses of nests to poaching
set.seed(20180515)
sim_len <- 400
N <- vector(mode = "list", length = sim_len)
nests = 100
p_poach <- 0.01
n <- c(50, 50, 50)
for (i in 1:sim_len) {
  if (nests > 0) {
    ## Probability of nest loss to poaching
    p <- runif(0, 1, n = nests)
    nests <- nests - length(which(p < p_poach))
  }
  A <- ysaFuncDD(yellow, n, nests)
  n <- A %*% n
  N[[i]] <- n
} 
map(N, sum) %>% unlist %>% plot

# 0.1% poaching 
## Simulation with losses of nests to poaching
set.seed(20180515)
sim_len <- 400
N <- vector(mode = "list", length = sim_len)
nests = 100
p_poach <- 0.1
n <- c(50, 50, 50)
for (i in 1:sim_len) {
  if (nests > 0) {
    ## Probability of nest loss to poaching
    p <- runif(0, 1, n = nests)
    nests <- nests - length(which(p < p_poach))
  }
  A <- ysaFuncDD(yellow, n, nests)
  n <- A %*% n
  N[[i]] <- n
}
map(N, sum) %>% unlist %>% plot

# 0.5% poaching 
## Simulation with losses of nests to poaching
set.seed(20180515)
sim_len <- 400
N <- vector(mode = "list", length = sim_len)
nests = 100
p_poach <- 0.5
n <- c(50, 50, 50)
for (i in 1:sim_len) {
  if (nests > 0) {
    ## Probability of nest loss to poaching
    p <- runif(0, 1, n = nests)
    nests <- nests - length(which(p < p_poach))
  }
  A <- ysaFuncDD(yellow, n, nests)
  n <- A %*% n
  N[[i]] <- n
}
map(N, sum) %>% unlist %>% plot

# 1% poaching 
## Simulation with losses of nests to poaching
set.seed(20180515)
sim_len <- 400
N <- vector(mode = "list", length = sim_len)
nests = 100
p_poach <- 1
n <- c(50, 50, 50)
for (i in 1:sim_len) {
  if (nests > 0) {
    ## Probability of nest loss to poaching
    p <- runif(0, 1, n = nests)
    nests <- nests - length(which(p < p_poach))
  }
  A <- ysaFuncDD(yellow, n, nests)
  n <- A %*% n
  N[[i]] <- n
}
map(N, sum) %>% unlist %>% plot

#So you can produce a plot which shows the abundance over time under different poaching scenarios. 
#And then a table or graph with the pop size at 50 and 100 years as a proportion of that Karrying Kapacity! 




