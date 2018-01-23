# Wisdom et al. 2000 
# life stage simulation analysis: estimating vital-rate effects on population growth for conservation 
library(dplyr)
library(ggplot2) 
library(knitr)
library(primer)
library(Rmisc)
library(agricolae) 
chicken <- read.csv("C:/Users/Ellie/OneDrive/Documents/1 UNIVERSITY/Level 4/Project & Dissertation/Wisdom 2000/chicken means sd .csv")
View(chicken)
tortoise <- read.csv("C:/Users/Ellie/OneDrive/Documents/1 UNIVERSITY/Level 4/Project & Dissertation/Wisdom 2000/desert tortoise means sd .csv")
View(tortoise)

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
tort_matrix

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

############ Application 1: simulated elasticity and regression metrics 
#1 random selection of vital rates 
#tortoise 
rnorm(100, mean = sx, sd = ssd) #x = mean, SD = standard deviation 

