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
library(tidyverse)
YSA_demog_data_master <- read.csv("YSA_demog_data_master.csv")
ysa <- YSA_demog_data_master 

View(ysa)

stage <- c("1a", "1b", "1c", "2", "3")
class <- c("egg", "nestling", "fledgling", "juvenile", "adult")
di <- c((27/356), (59/356), (270/356), (23/12), 7) #27 days, 59 days, To age 12 months, Age 13-36 months, Age 37 months+ (as 7 years)
pi <- c(0.89, 0.78, 0.73, 0.925, 0.925) #from meghann in YAS demography data csv 
F <- c(0, 0, 0, 0, 0.33)

yellow <- data_frame(stage, class, di, pi, F)

ysaFunc <- function (dataSource) 
{ 
#ps
p1a<- betaval((0.89), (0.06), fx=runif(1)) 
p1b<- betaval((0.78), (0.07), fx=runif(1))
p1c<- betaval((0.73), (0.07), fx=runif(1)) #0.07 as filler
p2 <- betaval((0.925), (0.025), fx=runif(1)) # 0.05 or 0.025?
p3 <- betaval((0.925), (0.025), fx=runif(1))
#f
f3 <- rnorm(1, mean = (3.2), sd = (0.24)) #should 3.3 be divided by 2  
#g 
Pi <- ((1 - (pi^(di - 1)))/(1 - (pi^di)))*pi #to include or not?
Gi <- (pi^di*(1 - pi))/(1 - pi^di)           #to include or not?
matrix2 <- matrix(0, nrow = 3, ncol = 3)
#add ps 
matrix2[1,1] <- (p1a*p1b*p1c)# this stage as the survival is from the multiplication of  p1a, p1b and p1c
matrix2[2,2] <- p2
matrix2[3,3] <- p3
#add f
matrix2[1,3] <- (f3)
#add gs 
matrix2[2,1] <- 
matrix2[3,2] <- Gi[4]
return(matrix2)
} 

ysaFunc()

## Example idea ---------------------------------

# make 10 matrices of random numbers.#for current code these aren't random and are drawn from ysa beta and normal distributed vital rates
mat1 <- map(1:10, ysaFunc)
# use these matrices in mat 1 and get the eigen system for each....
mat2<-map(mat1, function(x) eigen(x))

## Might work with yours... but needs your effort ---------------------------------

yellowMats <- map(1:10, function(x) yellowFunc(x)) #makes 10 matrices
outMats <- map(yellowMats, function(x) ysaFunc(x)) # makes 10 matrices
eigenOuts <- map(outMats, function(x) eigen(x))


############ for now 

yellowFunc <- function (dataSource)
{Pi <- ((1 - (pi^(di - 1)))/(1 - (pi^di)))*pi
Gi <- (pi^di*(1 - pi))/(1 - pi^di)
p1a <- Pi[1]
p1b <- Pi[2]
p1c <- Pi[3]
p2 <- Pi[4] 
p3 <- Pi[5]
f3 <- F[5]
p1 <- (pi[1]*pi[2]*pi[3])
d1 <- di[1]+di[2]+di[3] #because of 1a 1b 1c thing
P1 <- ((1 - (p1^(d1 - 1)))/(1 - (p1^d1)))*p1
g1 <- (p1^d1*(1 - p1))/(1 - p1^d1) #because of 1a 1b 1c thing
g2 <- Gi[4]
#matrix structure
matrix2 <- matrix(0, nrow = 3, ncol = 3)
#add ps 
matrix2[1,1] <- P1# this stage as the survival is from the multiplication of  p1a, p1b and p1c
matrix2[2,2] <- p2
matrix2[3,3] <- p3
#add f
matrix2[1,3] <- (f3)
#add gs 
matrix2[2,1] <- g1
matrix2[3,2] <- g2
return(matrix2)}

A <- yellowFunc(yellow)
A


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

#alter yellow so only have 3 stages 


#decresing fecundity and survival by 50%
yellow_fecAdjust<-mutate(yellow, F = 0.5*F)
yellow_pi1aAdj <- mutate(yellow, pi = ifelse(stage == "1a", pi * 0.5, pi * 1)) #funny (0.999 because it doesn't like 1s)
yellow_pi1bAdj <- mutate(yellow, pi = ifelse(stage == "1b", pi * 0.5, pi * 1))
yellow_pi1cAdj <- mutate(yellow, pi = ifelse(stage == "1c", pi * 0.5, pi * 1))
yellow_pi2Adj <- mutate(yellow, pi = ifelse(stage == "2", pi * 0.5, pi * 1))
yellow_pi3Adj <- mutate(yellow, pi = ifelse(stage == "3", pi * 0.5, pi * 1))
#50% increase in fecundity or an increase in survivorship to 1.0.
yellow_fecIncrease<-mutate(yellow, F = 1.5*F)
yellow_pi1aInc <- mutate(yellow, pi = ifelse(stage == "1a", (0.999-pi+pi), pi *1)) #funny (0.999 because it doesn't like 1s)
yellow_pi1bInc <- mutate(yellow, pi = ifelse(stage == "1b", (0.999-pi+pi), pi *1))
yellow_pi1cInc <- mutate(yellow, pi = ifelse(stage == "1c", (0.999-pi+pi), pi *1))
yellow_pi2Inc <- mutate(yellow, pi = ifelse(stage == "2", (0.999-pi+pi), pi *1))
yellow_pi3Inc <- mutate(yellow, pi = ifelse(stage == "3", (0.999-pi+pi), pi *1))
#decreases 
AFd <- yellowFunc(yellow_fecAdjust) 
eigs.AFd <- eigen(AFd)
A1ad <- yellowFunc(yellow_pi1aAdj) 
eigs.A1ad <- eigen(A1ad)
A1bd <- yellowFunc(yellow_pi1bAdj) 
eigs.A1bd <- eigen(A1bd)
A1cd <- yellowFunc(yellow_pi1cAdj) 
eigs.A1cd <- eigen(A1cd)
A2d <- yellowFunc(yellow_pi2Adj) 
eigs.A2d <- eigen(A2d)
A3d <- yellowFunc(yellow_pi3Adj) 
eigs.A3d <- eigen(A3d)
#increases 
AFi <- yellowFunc(yellow_fecIncrease) 
eigs.AFi <- eigen(AFi)
A1ai <- yellowFunc(yellow_pi1aInc) 
eigs.A1ai <- eigen(A1ai)
A1bi <- yellowFunc(yellow_pi1bInc) 
eigs.A1bi <- eigen(A1bi)
A1ci <- yellowFunc(yellow_pi1cInc) 
eigs.A1ci <- eigen(A1ci)
A2i <- yellowFunc(yellow_pi2Inc) 
eigs.A2i <- eigen(A2i)
A3i <- yellowFunc(yellow_pi3Inc) 
eigs.A3i <- eigen(A3i)
#finding the first eigenvalue (finite rate of increase)
#decreases
dom.pos.Fd <- which.max(eigs.AFd[["values"]])
LFd <- Re(eigs.AFd[["values"]][dom.pos.Fd]) 
dom.pos.1ad <- which.max(eigs.A1ad[["values"]])
L1ad <- Re(eigs.A1ad[["values"]][dom.pos.1ad]) 
dom.pos.1bd <- which.max(eigs.A1bd[["values"]])
L1bd <- Re(eigs.A1bd[["values"]][dom.pos.1bd]) 
dom.pos.1cd <- which.max(eigs.A1cd[["values"]])
L1cd <- Re(eigs.A1cd[["values"]][dom.pos.1cd]) 
dom.pos.2d <- which.max(eigs.A2d[["values"]])
L2d <- Re(eigs.A2d[["values"]][dom.pos.2d]) 
dom.pos.3d <- which.max(eigs.A3d[["values"]])
L3d <- Re(eigs.A3d[["values"]][dom.pos.3d]) 
#increases 
dom.pos.Fi <- which.max(eigs.AFi[["values"]])
LFi <- Re(eigs.AFi[["values"]][dom.pos.Fd]) 
dom.pos.1ai <- which.max(eigs.A1ai[["values"]])
L1ai <- Re(eigs.A1ai[["values"]][dom.pos.1ai]) 
dom.pos.1bi <- which.max(eigs.A1bi[["values"]])
L1bi <- Re(eigs.A1bi[["values"]][dom.pos.1bi]) 
dom.pos.1ci <- which.max(eigs.A1ci[["values"]])
L1ci <- Re(eigs.A1ci[["values"]][dom.pos.1ci]) 
dom.pos.2i <- which.max(eigs.A1ci[["values"]])
L2i <- Re(eigs.A2i[["values"]][dom.pos.2i]) 
dom.pos.3i <- which.max(eigs.A3i[["values"]])
L3i <- Re(eigs.A3i[["values"]][dom.pos.3i]) 

#plotting figure 1
lambdas<- c(LFd, L1ad, L1bd, L1cd, L2d, L3d)
rs <- log(lambdas)
class <- c("fecundity", "egg", "nestling", "fledgling", "juvenile", "adult")
table_decrease <- data.frame(class, rs)
graph <- ggplot(table_decrease, aes(x = class, y = rs)) + geom_bar(stat = "identity")
graph1 <- graph + labs(x = "Stage Class", y = "Intrinsic rate of Increase (r)") + linetype = "dashed", size = 0.75)  
graph2 <- graph1 + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ expand_limits(y = -0.25, x = 0.2)
figure1a <- graph2 + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + scale_x_discrete(limits = c("Fecundity", "Eggs/Hatchlings", "Small Juveniles", "Large Juveniles", "Subadults", "Novice Breeders", "1st-yr Remigrants", "Mature Breeders"))
figure1a 
class <- c("fecundity", "egg", "nestling", "fledgling", "juvenile", "adult")
lambdas<- c(LFi, L1ai, L1bi, L1ci, L2i, L3i)
rs <- log(lambdas)
table_decrease <- data.frame(class, rs)
graph <- ggplot(table_decrease, aes(x = class, y = rs)) + geom_bar(stat = "identity")
graph1 <- graph + labs(x = "Stage Class", y = "Intrinsic rate of Increase (r)") + geom_hline(aes(yintercept=-0.056399), linetype = "dashed", size = 0.75)  
graph2 <- graph1 + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ expand_limits(y = -0.25, x = 0.2)
figure1b <- graph2 + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + scale_x_discrete(limits = c("Fecundity", "Eggs/Hatchlings", "Small Juveniles", "Large Juveniles", "Subadults", "Novice Breeders", "1st-yr Remigrants", "Mature Breeders"))
figure1b 

