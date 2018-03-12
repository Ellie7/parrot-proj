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