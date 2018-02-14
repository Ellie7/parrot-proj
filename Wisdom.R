# Wisdom et al. 2000 
# life stage simulation analysis: estimating vital-rate effects on population growth for conservation 
rm(list=ls())
library(dplyr)
library(ggplot2) 
library(knitr)
library(primer)
library(Rmisc)
library(agricolae)
library(popbio)
library(MASS)
library(PerformanceAnalytics) #from online tutorial for rbeta
library(ExtDist)#from online tutorial #package ‘ExtDist’ is not available (for R version 3.4.2)
chicken <- read.csv("chicken means sd .csv")
View(chicken)
tortoise <- read.csv("desert tortoise means sd .csv")
View(tortoise)

################################################################ Steps 1-2 (as labelled in original paper)
# Application 1: simulated elasticity and regression metrics 
# making a function which produces a matrix based on randomly selected vital rates 
chickFunc <- function (chicken) 
{ sc1 <- betaval((chicken$s[1]), 0.129, fx=runif(1)) #0.129 is mean standard deviation for stage 1, didn't know how to correctly calculate ssd for
sc2 <- betaval((chicken$s[5]), (chicken$ssd[5]), fx=runif(1)) # this stage as the survival is from the multiplication of  s1a, s1b and s1c
sc3 <- betaval((chicken$s[6]), (chicken$ssd[6]), fx=runif(1))
sc4 <- betaval((chicken$s[7]), (chicken$ssd[7]), fx=runif(1))
#m 
mc2 <- rnorm(1, mean = (chicken$m[5]), sd = (chicken$msd[5]))
mc3 <- rnorm(1, mean = (chicken$m[6]), sd = (chicken$msd[6]))
mc4 <- rnorm(1, mean = (chicken$m[7]), sd = (chicken$msd[7]))
mc5 <- rnorm(1, mean = (chicken$m[7]), sd = (chicken$msd[7]))
matrix2 <- matrix(0, nrow = 4, ncol = 4)
#add sxmx
matrix2[1,1] <- (sc1*mc2) 
matrix2[1,2] <- (sc2*mc3)
matrix2[1,3] <- (sc3*mc4)
matrix2[1,4] <- (sc4*mc5)
#add s(mxc$m[7]))
matrix2[2,1] <- sc1
matrix2[3,2] <- sc2
matrix2[4,3] <- sc3
matrix2[4,4] <- sc4
return(matrix2)
}

chickFunc(chicken) #chickFunc returns a different matrix each time 

tortFunc <- function(tortoise)
{ #s
  s1 <- betaval((tortoise$s[1]), (tortoise$ssd[1]), fx=runif(1)) # betaval returns a random beta value 
  s2 <- betaval((tortoise$s[2]), (tortoise$ssd[2]), fx=runif(1)) # usage: betaval(mn, sdev, fx=runif(1))
  s3 <- betaval((tortoise$s[3]), (tortoise$ssd[3]), fx=runif(1))
  s4 <- betaval((tortoise$s[4]), (tortoise$ssd[4]), fx=runif(1))
  s5 <- betaval((tortoise$s[5]), (tortoise$ssd[5]), fx=runif(1))
  s6 <- betaval((tortoise$s[6]), (tortoise$ssd[6]), fx=runif(1))
  s7 <- betaval((tortoise$s[7]), (tortoise$ssd[7]), fx=runif(1))
  s8 <- betaval((tortoise$s[8]), (tortoise$ssd[8]), fx=runif(1))
  #g 
  g2 <- betaval((tortoise$g[2]), (tortoise$gsd[2]), fx=runif(1))
  g3 <- betaval((tortoise$g[3]), (tortoise$gsd[3]), fx=runif(1))
  g4 <- betaval((tortoise$g[4]), (tortoise$gsd[4]), fx=runif(1))
  g5 <- betaval((tortoise$g[5]), (tortoise$gsd[5]), fx=runif(1))
  g6 <- betaval((tortoise$g[6]), (tortoise$gsd[6]), fx=runif(1))
  g7 <- betaval((tortoise$g[7]), ((tortoise$gsd[7])/2), fx=runif(1)) # /2 because SD too large for beta distribution 
  #m
  m6 <- rnorm(1, mean = (tortoise$m[6]), sd = (tortoise$msd[6]))
  m7 <- rnorm(1, mean = (tortoise$m[7]), sd = (tortoise$msd[7]))
  m8 <- rnorm(1, mean = (tortoise$m[8]), sd = (tortoise$msd[8]))
  #make matrix
  matrix1 <- matrix(0, nrow = 8, ncol = 8)
  #add Fs
  matrix1[1,6] <- m6
  matrix1[1,7] <- m7
  matrix1[1,8] <- m8
  #add Ps (diagonals)
  matrix1[2,2] <- (s2*(1-g2)) 
  matrix1[3,3] <- (s3*(1-g3)) 
  matrix1[4,4] <- (s4*(1-g4)) 
  matrix1[5,5] <- (s5*(1-g5)) 
  matrix1[6,6] <- (s6*(1-g6)) 
  matrix1[7,7] <- (s7*(1-g7)) 
  matrix1[8,8] <- s8 
  matrix1
  #add Gs (off-diagonals)
  matrix1[2,1] <- s1
  matrix1[3,2] <- s2*g2
  matrix1[4,3] <- s3*g3
  matrix1[5,4] <- s4*g4
  matrix1[6,5] <- s5*g5
  matrix1[7,6] <- s6*g6
  matrix1[8,7] <- s7*g7
  #remove NAs
  #shows location of NAs
  is.na(matrix1)
  #replaces NAs with 0s
  matrix1[is.na(matrix1)] <- 0
  matrix1#almost there 
  return(matrix1)}

tortFunc(tortoise) #tortFunc returns a different matrix each time drawn from beta & lognorm distributions 
#currently erroring: "Error during wrapup: Standard deviation too high for beta distribution" 

# Generating correlated vital rates using an estimated correlation between matrix between vital rates 
# using chapter 8, box 8.6 matlab code from Morris & Doak 2002 
#----- simulation parameters
#parameters for two vital rates  (s and m)
# a beta and a lognormal 
vrmeans <- c(0.0945, 0.445, 0.51, 0.284) # means for survival only atm, in book method is for three vital rates 
vrvars <- c(0, (0.081^2), (0.079^2), (0.090^2)) # variances #^2 because standard deviation is the square root of the variance 

#minimum and maximum values for each vital rate 
# zeros are palceholders for rates that are not stretched betas
vrmins <- c(0, 0, 0, 0)
vrmaxs <- c(0, 0, 0, 0) 
#then a full correlation matrix 
cor(vrmeans, vrvars)
#find the number of vital rates

#find the eigen values (D) and eigenvectors (W) of the correlation matrix 

#calculate C12 the marix to use to make correlated standard normal variates from uncorrelated ones
C12 <- W*(sqrt(abs(D)))*W 
#loop to do each years vita rates 

mvrnorm(n = 1, mu, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE) # from online 

########################################################################### Step 3
#the process is excecuted 1000 times, resulting in 1000 matrix replicates
# of vital rates and matrix elements for a given species 



############################################################################ Step 4 
#lambda and lower level elasticities associated with each vital rate 
# calculated at a stable stage ditribution

A <- matrix
#eigen analysis 
eigs.A <- eigen(A)
eigs.A
#finding the first eigenvalue (finite rate of increase)
dom.pos <- which.max(eigs.A[["values"]])
L1 <- Re(eigs.A[["values"]][dom.pos])
L1
lambda <- Re(eigs.A$values[1])

#finding r 
r <- log(L1)
r

#calculating the stable stage distribution 
w <- Re(eigs.A[["vectors"]][, dom.pos])
ssd <- w / sum(w)
stable <- ssd*100

############################################################################ Step 5
#Data across replicates were analyzed to estimate effects of each 
# vital rate on lambda 

