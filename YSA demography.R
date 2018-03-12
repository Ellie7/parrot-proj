# Yellow-shouldered amazon  
rm(list=ls())
library(dplyr)
library(ggplot2) 
library(knitr)
library(primer)
library(Rmisc)
library(agricolae)
library(popbio)
library(MASS) 
YSA_demog_data_master <- read.csv("YSA_demog_data_master.csv")
ysa <- YSA_demog_data_master 



ysaFunc <- function (ysa) 
{ #ps
  p1 <- betaval((chicken$s[1]), 0.129, fx=runif(1)) #0.129 is mean standard deviation for stage 1, didn't know how to correctly calculate ssd for
p2 <- betaval((chicken$s[5]), (chicken$ssd[5]), fx=runif(1)) # this stage as the survival is from the multiplication of  s1a, s1b and s1c
p3 <- betaval((chicken$s[6]), (chicken$ssd[6]), fx=runif(1))
#f
f3 <- rnorm(1, mean = (chicken$m[7]), sd = (chicken$msd[7]))
#g 
g1 <- 
g2 <- 
matrix2 <- matrix(0, nrow = 4, ncol = 4)
#add first row
matrix2[1,1] <- (p1a*p1b*p1c) 
matrix2[1,2] <- (f2)
matrix2[1,3] <- (f3)
#add p
matrix2[2,1] <- g1
matrix2[3,2] <- p2
matrix2[4,3] <- p3
return(matrix2)
}