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
  

{ 
#ps
p1a<- betaval((chicken$s[1]), 0.129, fx=runif(1)) 
 
p2 <- betaval((chicken$s[5]), (chicken$ssd[5]), fx=runif(1)) 
p3 <- betaval((chicken$s[6]), (chicken$ssd[6]), fx=runif(1))
#f
f3 <- rnorm(1, mean = (chicken$m[7]), sd = (chicken$msd[7]))
#g 
g1 <- 
g2 <- 
matrix2 <- matrix(0, nrow = 4, ncol = 4)
#add ps 
matrix2[1,1] <- (p1a*p1b*p1c)# this stage as the survival is from the multiplication of  p1a, p1b and p1c
matrix2[2,2] <- p2
matrix2[3,3] <- p3
#add f
matrix2[1,3] <- (f3)
#add gs 
matrix2[2,1] <- g1
matrix2[3,2] <- g2
return(matrix2)
}