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
#finding out carrying capacity for 100 nests 
sim_len <- 200
N <- vector(mode = "list", length = sim_len)
thres = 100
n <- c(50, 50, 50)
for (i in 1:sim_len) {
  A <- ysaFuncDD(yellow, n, thres)
  n <- A %*% n
  N[[i]] <- n
}
map(N, sum) %>% unlist %>% plot

number <- map(N, sum) %>% unlist
index <- c(1:200) 

K100 <- max(number)

#finding out carrying capacity for 200 nests 
sim_len <- 200
N <- vector(mode = "list", length = sim_len)
thres = 200
n <- c(50, 50, 50)
for (i in 1:sim_len) {
  A <- ysaFuncDD(yellow, n, thres)
  n <- A %*% n
  N[[i]] <- n
}
map(N, sum) %>% unlist %>% plot

number2 <- map(N, sum) %>% unlist
index <- c(1:200) 

K200 <- max(number2)

dens <- data_frame(number, index)

dens_plot <- ggplot(dens, aes(x = index, y = number)) + geom_point() + mytheme +
  labs(x = "Time", y = "Population Size") 
dens_plot + xlim(0,100)

#carrying capacity for 100 nests 
K100
#carrying capacity for 200 nests 
K200 
#stable stage disribution
S <- ysaFunc(yellow)
eigen_S <- eigen.analysis(S)
stable_stage <- c(eigen_S$stable.stage)
stable_stage <- data_frame(stable_stage)
stable_stage
#stable stage for 100 nests 

stable_100 <- c((K100 * stable_stage$stable_stage[1]),(K100 * stable_stage$stable_stage[2]), (K100 * stable_stage$stable_stage[3]))
#stable stage for 200 nests
stable_200 <- c((K200 * stable_stage$stable_stage[1]),(K200 * stable_stage$stable_stage[2]), (K200 * stable_stage$stable_stage[3]))

#-----------------------------------------------------------------------------------------------------------------------------------
#use an implementation of the function described above to simulate population trajectories under a range of poaching 
# probabilities: .01%, 0.1%, 0.5%, 1%
# used to simulate the population size under different levels of poaching related nest site loss at say 50 or 100 years in the 
#future.
# 0.01% poaching 
## Simulation with losses of nests to poaching 
set.seed(20180515)
sim_len <- 100
Na <- vector(mode = "list", length = sim_len)
nests = 100
p_poach <- 0.0001
n <- stable_100 
for (i in 1:sim_len) {
  if (nests > 0) {
    ## Probability of nest loss to poaching
    p <- runif(0, 1, n = nests)
    nests <- nests - length(which(p < p_poach))
  }
  A <- ysaFuncDD(yellow, n, nests)
  n <- A %*% n
  Na[[i]] <- n
} 
map(Na, sum) %>% unlist %>% plot

number <- map(Na, sum) %>% unlist
index <- c(1:100) 

dens <- data_frame(number, index)

dens_plot <- ggplot(dens, aes(x = index, y = number)) + geom_point() + mytheme +
  labs(x = "Time", y = "Population Size") 
dens_plot + xlim(0,100)

# 0.1% poaching 
## Simulation with losses of nests to poaching
set.seed(20180515)
sim_len <- 100
Nb <- vector(mode = "list", length = sim_len)
nests = 100
p_poach <- 0.001
n <- stable_100
for (i in 1:sim_len) {
  if (nests > 0) {
    ## Probability of nest loss to poaching
    p <- runif(0, 1, n = nests)
    nests <- nests - length(which(p < p_poach))
  }
  A <- ysaFuncDD(yellow, n, nests)
  n <- A %*% n
  Nb[[i]] <- n
}
map(Nb, sum) %>% unlist %>% plot

# 0.5% poaching 
## Simulation with losses of nests to poaching
set.seed(20180515)
sim_len <- 100
Nc <- vector(mode = "list", length = sim_len)
nests = 100
p_poach <- 0.005
n <- stable_100
for (i in 1:sim_len) {
  if (nests > 0) {
    ## Probability of nest loss to poaching
    p <- runif(0, 1, n = nests)
    nests <- nests - length(which(p < p_poach))
  }
  A <- ysaFuncDD(yellow, n, nests)
  n <- A %*% n
  Nc[[i]] <- n
}
map(Nc, sum) %>% unlist %>% plot

# 1% poaching 
## Simulation with losses of nests to poaching
set.seed(20180515)
sim_len <- 100
Nd <- vector(mode = "list", length = sim_len)
nests = 100
p_poach <- 0.01
n <- stable_100
for (i in 1:sim_len) {
  if (nests > 0) {
    ## Probability of nest loss to poaching
    p <- runif(0, 1, n = nests)
    nests <- nests - length(which(p < p_poach))
  }
  A <- ysaFuncDD(yellow, n, nests)
  n <- A %*% n
  Nd[[i]] <- n
}
map(Nd, sum) %>% unlist %>% plot

#So you can produce a plot which shows the abundance over time under different poaching scenarios. 
#And then a table or graph with the pop size at 50 and 100 years as a proportion of that Karrying Kapacity! 

Na <- map(Na, sum) %>% unlist
Nb <- map(Nb, sum) %>% unlist
Nc <- map(Nc, sum) %>% unlist
Nd <- map(Nd, sum) %>% unlist

Ns <- c(Na, Nb, Nc, Nd)
Ind <- c(1:100, 1:100, 1:100, 1:100)
Poaching <- c(rep("0.01%", 100),rep("0.1%", 100), rep("1.0%", 100),rep("5.0%", 100)) 
poach <- data_frame(Ns, Ind, Poaching)

poach_plot <- ggplot(poach, aes(x = Ind, y = Ns, group_by(Poaching))) + geom_point(aes(colour = Poaching)) + mytheme +
  labs(x = "Time", y = "Population Size") 
poach_plot

#-----------------------------------------------------------------------------------------------------------------------------------
# if there are 200 nests 
# 0.01% poaching 
## Simulation with losses of nests to poaching
set.seed(20180515)
sim_len <- 100
Na2 <- vector(mode = "list", length = sim_len)
nests = 200
p_poach <- 0.0001
n <- stable_200
for (i in 1:sim_len) {
  if (nests > 0) {
    ## Probability of nest loss to poaching
    p <- runif(0, 1, n = nests)
    nests <- nests - length(which(p < p_poach))
  }
  A <- ysaFuncDD(yellow, n, nests)
  n <- A %*% n
  Na2[[i]] <- n
} 
map(Na2, sum) %>% unlist %>% plot

# 0.1% poaching 
## Simulation with losses of nests to poaching
set.seed(20180515)
sim_len <- 100
Nb2 <- vector(mode = "list", length = sim_len)
nests = 200
p_poach <- 0.001
n <- stable_200
for (i in 1:sim_len) {
  if (nests > 0) {
    ## Probability of nest loss to poaching
    p <- runif(0, 1, n = nests)
    nests <- nests - length(which(p < p_poach))
  }
  A <- ysaFuncDD(yellow, n, nests)
  n <- A %*% n
  Nb2[[i]] <- n
}
map(Nb2, sum) %>% unlist %>% plot

# 0.5% poaching 
## Simulation with losses of nests to poaching
set.seed(20180515)
sim_len <- 100
Nc2 <- vector(mode = "list", length = sim_len)
nests = 200
p_poach <- 0.005
n <- stable_200
for (i in 1:sim_len) {
  if (nests > 0) {
    ## Probability of nest loss to poaching
    p <- runif(0, 1, n = nests)
    nests <- nests - length(which(p < p_poach))
  }
  A <- ysaFuncDD(yellow, n, nests)
  n <- A %*% n
  Nc2[[i]] <- n
}
map(Nc2, sum) %>% unlist %>% plot

# 1% poaching 
## Simulation with losses of nests to poaching
set.seed(20180515)
sim_len <- 100
Nd2 <- vector(mode = "list", length = sim_len)
nests = 200
p_poach <- 0.01
n <- stable_200
for (i in 1:sim_len) {
  if (nests > 0) {
    ## Probability of nest loss to poaching
    p <- runif(0, 1, n = nests)
    nests <- nests - length(which(p < p_poach))
  }
  A <- ysaFuncDD(yellow, n, nests)
  n <- A %*% n
  Nd2[[i]] <- n
}
map(Nd2, sum) %>% unlist %>% plot

#So you can produce a plot which shows the abundance over time under different poaching scenarios. 
#And then a table or graph with the pop size at 50 and 100 years as a proportion of that Karrying Kapacity! 

Na2 <- map(Na2, sum) %>% unlist
Nb2 <- map(Nb2, sum) %>% unlist
Nc2 <- map(Nc2, sum) %>% unlist
Nd2 <- map(Nd2, sum) %>% unlist

Ns2 <- c(Na2, Nb2, Nc2, Nd2)
Ind <- c(1:100, 1:100, 1:100, 1:100)
Poaching <- c(rep("0.01%", 100),rep("0.1%", 100), rep("1.0%", 100),rep("5.0%", 100)) 
poach2 <- data_frame(Ns2, Ind, Poaching)

poach_plot_200 <- ggplot(poach, aes(x = Ind, y = Ns2, group_by(Poaching))) + geom_point(aes(colour = Poaching)) + mytheme +
  labs(x = "Time", y = "Population Size") 
poach_plot_200

#-----------------------------------------------------------------------------------------------------------------------------------
#facet 
Nsf <- c(Na, Nb, Nc, Nd,Na2, Nb2, Nc2, Nd2)
Indf <- c(1:100, 1:100, 1:100, 1:100, 1:100, 1:100, 1:100, 1:100)
Poaching <- c(rep("0.01%", 100),rep("0.1%", 100), rep("0.5%", 100),rep("1.0%", 100),
              rep("0.01%", 100),rep("0.1%", 100), rep("0.5%", 100),rep("1.0%", 100)) 
Nests_av <- c(rep("100 Nest Sites",400), rep("200 Nest Sites", 400))
stable_pop <- c(rep(K100, 400), rep(K200, 400))
poach_facet_df <- data_frame(Nsf, Indf, Poaching, Nests_av, stable_pop)

poach_facet <- ggplot(poach_facet_df, aes(x = Indf, y = Nsf, group_by(Poaching))) + geom_point(aes(colour = Poaching)) + facet_wrap(~ Nests_av)
poach_facet<- poach_facet + mytheme + labs(x = "Time (years)", y = "Population Size") 
poach_facet + geom_vline(xintercept = 50, linetype = "dashed")

#-----------------------------------------------------------------------------------------------------------------------------------
# make a table (gasp) of the %decrease from stable pop at 50 or 100 years for each sim.
poach_facet_df

gasp <- subset(poach_facet_df, Indf == 50, select = c(Nsf, Poaching, Nests_av, stable_pop))

gasp <- mutate(gasp, decrease = (((stable_pop-Nsf)/stable_pop)*100))

tab_nest <- subset(gasp, select = c(Poaching, Nests_av, decrease))

colnames(tab_nest) <- c("Poaching", "Threshold", "% Decrease from stable population size 50 years on" )
tab_nest
kable(tab_nest, caption = "Table x. Percentage decrease in population size 50 years into the future for two different population 
      thresholds, 100 nest sites and 200 nest sites.")
