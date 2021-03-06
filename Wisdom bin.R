# Wisdom bin 
library(dplyr)
library(ggplot2) 
library(knitr)
library(primer)
library(Rmisc)
library(agricolae)
library(popbio)
library(PerformanceAnalytics) #from online tutorial 
library(ExtDist)#from online tutorial #package ‘ExtDist’ is not available (for R version 3.4.2)
chicken <- read.csv("chicken means sd .csv")
View(chicken)
tortoise <- read.csv("desert tortoise means sd .csv")
View(tortoise)
---------------------------------------------------------------- #wisdom bin 
#Creating a stage-based projection matrix for desert tortoise 
lifetable <- tortoise
tortoiseFunc <- function (lifetable) 
{
  sx <- select(lifetable, s)
  gx <- select(lifetable, g)
  diags <- (sx*(1-gx))
  offdiags <- (sx*gx)
  matrix1 <- matrix(0, nrow = 8, ncol = 8)
  #add Fs
  matrix1[1,] <- lifetable$m 
  #add Ps (diagonals)
  matrix1[2,2] <- diags$s[2] 
  matrix1[3,3] <- diags$s[3]
  matrix1[4,4] <- diags$s[4]
  matrix1[5,5] <- diags$s[5]
  matrix1[6,6] <- diags$s[6]
  matrix1[7,7] <- diags$s[7]
  matrix1[8,8] <- diags$s[8]
  matrix1
  #add Gs (off-diagonals)
  matrix1[2,1] <- sx$s[1]
  matrix1[3,2] <- offdiags$s[2]
  matrix1[4,3] <- offdiags$s[3]
  matrix1[5,4] <- offdiags$s[4]
  matrix1[6,5] <- offdiags$s[5]
  matrix1[7,6] <- offdiags$s[6] 
  matrix1[8,7] <- offdiags$s[7] 
  #remove NAs
  #shows location of NAs
  is.na(matrix1)
  #replaces NAs with 0s
  matrix1[is.na(matrix1)] <- 0
  matrix1#almost there 
  return(matrix1)
}

tort_matrix <- tortoiseFunc(tortoise) 
tort_matrix #matrix with mean valyes 

#Creating a stage-based projection matrix for prairie chicken 
lifetable <- chicken 
chickenFunc <- function (lifetable) 
{
  sxc <- select(lifetable, s) #c in sxc and mxc stands for chicken 
  mxc <- select(lifetable, m)
  matrix2 <- matrix(0, nrow = 4, ncol = 4)
  #add sxmx
  matrix2[1,1] <- ((sxc$s[1])*(mxc$m[5])) #numbers look weird because data is in 8 age classes, this matrix is 4x4 (see appendix pg 641)
  matrix2[1,2] <- ((sxc$s[5])*(mxc$m[6])) 
  matrix2[1,3] <- ((sxc$s[6])*(mxc$m[7]))
  matrix2[1,4] <- ((sxc$s[7])*(mxc$m[7]))
  #add s
  matrix2[2,1] <- sxc$s[1]
  matrix2[3,2] <- sxc$s[5] #numbers look weird because data is in 8 age classes, this matrix is 4x4 (see appendix pg 641)
  matrix2[4,3] <- sxc$s[6]
  matrix2[4,4] <- sxc$s[7]
  matrix2
  return(matrix2)
}

chickenFunc(chicken) 

#1 random selection of vital rates 
#tortoise 
#s
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
g7 <- betaval((tortoise$g[7]), (tortoise$gsd[7]), fx=runif(1))
#m
m6 <- rnorm(100, mean = (tortoise$m[6]), sd = (tortoise$msd[6]))
m7 <- rnorm(100, mean = (tortoise$m[7]), sd = (tortoise$msd[7]))
m8 <- rnorm(100, mean = (tortoise$m[8]), sd = (tortoise$msd[8]))

#chicken 
#s
sc1 <- betaval((chicken$s[1]), 0.129, fx=runif(1)) #0.129 is mean standard deviation for stage 1, didn't know how to correctly calculate ssd for
sc2 <- betaval((chicken$s[5]), (chicken$ssd[5]), fx=runif(1)) # this stage as the survival is from the multiplication of  s1a, s1b and s1c
sc3 <- betaval((chicken$s[6]), (chicken$ssd[6]), fx=runif(1))
sc4 <- betaval((chicken$s[7]), (chicken$ssd[7]), fx=runif(1))
#m 
mc2 <- rnorm(100, mean = (chicken$m[5]), sd = (chicken$msd[5]))
mc3 <- rnorm(100, mean = (chicken$m[6]), sd = (chicken$msd[6]))
mc4 <- rnorm(100, mean = (chicken$m[7]), sd = (chicken$msd[7]))
mc5 <- rnorm(100, mean = (chicken$m[7]), sd = (chicken$msd[7]))