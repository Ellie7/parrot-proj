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
di <- c((27/356), (59/356), (270/356), (23/12), 7) #27 days, 59 days, To age 12 months, Age 13-36 months, Age 37 months+ (as 7 years)
pi <- c(0.89, 0.78, 0.73, 0.925, 0.925) #from meghann in YAS demography data csv 
F <- c(0, 0, 0, 0, 0.33)

yellow <- data_frame(stage, class, di, pi, F)

ysaFunc <- function (ysa) 
{ 
#ps
p1a<- betaval(), (), fx=runif(1)) 
 
p2 <- betaval(), (), fx=runif(1)) 
p3 <- betaval(), (), fx=runif(1))
#f
f3 <- rnorm(1, mean = (0.33), sd = ()) #should 0.33 be divided by 2? because its the percentage of females 
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

############ for now 

yellowFunc <- function (yellow)
{Pi <- ((1 - (pi^(di - 1)))/(1 - (pi^di)))*pi
Gi <- (pi^di*(1 - pi))/(1 - pi^di)
p1a <- Pi[1]
p1b <- Pi[2]
p1c <- Pi[3]
p2 <- Pi[4] 
p3 <- Pi[5]
f3 <- F[5]
g1 <- Gi[1]
g2 <- Gi[2]
#matrix structure
matrix2 <- matrix(0, nrow = 3, ncol = 3)
#add ps 
matrix2[1,1] <- (p1a*p1b*p1c)# this stage as the survival is from the multiplication of  p1a, p1b and p1c
matrix2[2,2] <- p2
matrix2[3,3] <- p3
#add f
matrix2[1,3] <- (f3)
#add gs 
matrix2[2,1] <- g1
matrix2[3,2] <- g2
A <- matrix2
A}

A <- yellowFunc(yellow)


#eigen analysis 
eigs.A <- eigen(A)
eigs.A
#finding the first eigenvalue (finite rate of increase)
dom.pos <- which.max(eigs.A[["values"]])
L1 <- Re(eigs.A[["values"]][dom.pos])
L1
lambda <- Re(eigs.A$values[1])
#=0.9329156
#finding r 
r <- log(L1)
r
#=0.06944057
w <- Re(eigs.A[["vectors"]][, dom.pos])
ssd <- w / sum(w)
stable <- ssd*100
stable <- round(stable, 2)
#0.207 0.670 0.115 0.007 0.000 0.000 0.002 
#calculating the reproductive value 
M <- eigen(t(A))
v <- Re(M$vectors[,which.max(Re(M$values))])
RV <- v / v[1]
RV 
#creating table 
stage_number <- c(1,2,3)
class_label <- c("Egg", "Juvenile", "Adult")
tab_5<- data.frame(stage_number, class_label, stable, RV)
colnames(tab_5) <- c("Stage number", "Stage Class", "Stable stage distribution (Dominant eigenvector)", "Reproductive values (left eigenvector)" )
tab_5
kable(tab_5, caption = "Table 1. Stable stage distribution (wJ) and reproductive values (v') for the yellow shouldered amazon population matrix.")
