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
  matrix1[1,1] <- 0 
  matrix1[2,2] <- diags$s[2] 
  matrix1[3,3] <- diags$s[3]
  matrix1[4,4] <- diags$s[4]
  matrix1[5,5] <- diags$s[5]
  matrix1[6,6] <- diags$s[6]
  matrix1[7,7] <- diags$s[7]
  matrix1[8,8] <- diags$s[8]
  matrix1
  #add Gs (off-diagonals)
  matrix1[2,1] <- sx[1]
  matrix1[3,2] <- offdiags$s[2]
  matrix1[4,3] <- offdiags$s[3]
  matrix1[5,4] <- offdiags$s[4]
  matrix1[6,5] <- offdiags$s[5]
  matrix1[7,6] <- offdiags$s[6] 
  matrix1[8,7] <- offdiags$s[7] 
  return(matrix1)
}

tortoiseFunc(tortoise) 
