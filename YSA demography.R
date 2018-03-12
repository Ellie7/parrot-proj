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

View(ysa)

stage <- c("1a", "1b", "1c", "2", "3")
class <- c("egg", "nestling", "fledgling", "juvenile", "adult")
stage_duration <- ((27/356), (59/356), (270/356), (23/12), 7) #27 days, 59 days, To age 12 months, Age 13-36 months, Age 37 months+ (as 7 years)
 
p <- c()
g <- c()
ysaFunc <- function (ysa) 
{ 
#ps
p1a<- betaval(), (), fx=runif(1)) 
 
p2 <- betaval(), (), fx=runif(1)) 
p3 <- betaval(), (), fx=runif(1))
#f
f3 <- rnorm(1, mean = (), sd = ())
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