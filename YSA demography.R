# Yellow-shouldered amazon elasticity analyses 
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
source(file = "ysa functions.R")

YSA_demog_data_master <- read.csv("YSA_demog_data_master.csv")
ysa <- YSA_demog_data_master 

stage <- c("1a", "1b", "1c", "2", "3")
class <- c("egg", "nestling", "fledgling", "juvenile", "adult")
di <- c((27/356), (59/356), (270/356), (23/12), 7) #27 days, 59 days, To age 12 months, Age 13-36 months, Age 37 months+ (as 7 years)
pi <- c(0.89, 0.78, 0.73, 0.925, 0.925) #from meghann in YAS demography data csv 
piSD <- c(0.06, 0.07, 0.07, 0.025, 0.025) # 2nd 0.07 as filler for now 
f <- c(0, 0, 0, 0, 3.2) #may need to halve 3.2 (sex ratio assumed 1:1)
fSD <- c(0, 0, 0, 0, 0.24) # may need to halve 

yellow <- data_frame(stage, class, di, pi, piSD,  f, fSD)

## Example idea ---------------------------------
ysaFunc(yellow)
yellow
#make 10 matrices of random numbers.#for current code these aren't random and are drawn from ysa beta and normal distributed vital rates
mat1 <- map(1:10, function(x) ysaFunc(yellow))
# use these matrices in mat 1 and get the eigen system for each....
mat2<-map(mat1, function(x) eigen(x))

A <- ysameanFunc(yellow) # for 'mean' matrix 
A <- ysaFunc(yellow) #for matrix drawn randomly from beta and lognormal distributed vital rates 
#eigen analysis 
eigs.A <- eigen(A)
eigs.A
#finding the first eigenvalue (finite rate of increase)
dom.pos <- which.max(eigs.A[["values"]])
L1mean <- Re(eigs.A[["values"]][dom.pos])
L1mean
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

#decresing fecundity and survival by 10%
yellow_fecAdjust<-mutate(yellow, F = 0.9*F)
yellow_pi1aAdj <- mutate(yellow, pi = ifelse(stage == "1a", pi * 0.9, pi * 1)) 
yellow_pi1bAdj <- mutate(yellow, pi = ifelse(stage == "1b", pi * 0.9, pi * 1))
yellow_pi1cAdj <- mutate(yellow, pi = ifelse(stage == "1c", pi * 0.9, pi * 1))
yellow_pi2Adj <- mutate(yellow, pi = ifelse(stage == "2", pi * 0.9, pi * 1))
yellow_pi3Adj <- mutate(yellow, pi = ifelse(stage == "3", pi * 0.9, pi * 1))
#increase in fecundity or an increase in survivorship by 10%
yellow_fecIncrease<-mutate(yellow, F = 1.5*F)
yellow_pi1aInc <- mutate(yellow, pi = ifelse(stage == "1a", pi * 1.1, pi *1))
yellow_pi1bInc <- mutate(yellow, pi = ifelse(stage == "1b", pi * 1.1, piSD *1))
yellow_pi1cInc <- mutate(yellow, pi = ifelse(stage == "1c", pi * 1.1, pi *1))
yellow_pi2Inc <- mutate(yellow, pi = ifelse(stage == "2", pi * 1.05, pi *1))
yellow_pi3Inc <- mutate(yellow, pi = ifelse(stage == "3", pi * 1.05, pi *1))
#decreases 
AFd <- ysaFunc(yellow_fecAdjust) 
eigs.AFd <- eigen(AFd)
A1ad <- ysaFunc(yellow_pi1aAdj) 
eigs.A1ad <- eigen(A1ad)
A1bd <- ysaFunc(yellow_pi1bAdj) 
eigs.A1bd <- eigen(A1bd)
A1cd <- ysaFunc(yellow_pi1cAdj) 
eigs.A1cd <- eigen(A1cd)
A2d <- ysaFunc(yellow_pi2Adj) 
eigs.A2d <- eigen(A2d)
A3d <- ysaFunc(yellow_pi3Adj) 
eigs.A3d <- eigen(A3d)
#increases 
AFi <- ysaFunc(yellow_fecIncrease) 
eigs.AFi <- eigen(AFi)
A1ai <- ysaFunc(yellow_pi1aInc) 
eigs.A1ai <- eigen(A1ai)
A1bi <- ysaFunc(yellow_pi1bInc) 
eigs.A1bi <- eigen(A1bi)
A1ci <- ysaFunc(yellow_pi1cInc) 
eigs.A1ci <- eigen(A1ci)
A2i <- ysaFunc(yellow_pi2Inc) 
eigs.A2i <- eigen(A2i)
A3i <- ysaFunc(yellow_pi3Inc) 
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
LFi <- Re(eigs.AFi[["values"]][dom.pos.Fi]) 
dom.pos.1ai <- which.max(eigs.A1ai[["values"]])
L1ai <- Re(eigs.A1ai[["values"]][dom.pos.1ai]) 
dom.pos.1bi <- which.max(eigs.A1bi[["values"]])
L1bi <- Re(eigs.A1bi[["values"]][dom.pos.1bi]) 
dom.pos.1ci <- which.max(eigs.A1ci[["values"]])
L1ci <- Re(eigs.A1ci[["values"]][dom.pos.1ci]) 
dom.pos.2i <- which.max(eigs.A2i[["values"]])
L2i <- Re(eigs.A2i[["values"]][dom.pos.2i]) 
dom.pos.3i <- which.max(eigs.A3i[["values"]])
L3i <- Re(eigs.A3i[["values"]][dom.pos.3i]) 

#plotting a figure -  Changes in rate of increase r resulting from simulated changes (10%) in fecundity and survival of individual 
#life history stages  (similar to figure 1a & 1b in crouse 1987)
lambdas<- c(LFd, L1ad, L1bd, L1cd, L2d, L3d)
rsa <- log(lambdas)
class <- c("fecundity", "egg", "nestling", "fledgling", "juvenile", "adult")
table_decrease <- data.frame(class, rsa)
grapha <- ggplot(table_decrease, aes(x = class, y = rsa)) + geom_bar(stat = "identity")
graph1a <- grapha + labs(x = "Stage Class", y = "Intrinsic rate of Increase (r)") 
graph2a <- graph1a + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())
figure1a <- graph2a + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + scale_x_discrete(limits = class)
figure1a 
lambdas<- c(LFi, L1ai, L1bi, L1ci, L2i, L3i)
rsb <- log(lambdas)
table_decrease <- data.frame(class, rsb)
graphb <- ggplot(table_decrease, aes(x = class, y = rsb)) + geom_bar(stat = "identity")
graph1b <- graphb + labs(x = "Stage Class", y = "Intrinsic rate of Increase (r)") 
graph2b <- graph1b + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())
figure1b <- graph2b + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + scale_x_discrete(limits = class)
figure1b 

######## figure 3 
#the elasticity, or proportional sensitivity of lambda to changes in fecundity F, survival while remaining in the same stage P, and survival 
# with growth, G. Because the elasticities of these matrix elements sum to 1, they can be compared directly in terms of their contribution to the 
# population growth rate 
A <- ysameanFunc(yellow)
eigs.A <- eigen(A)
eigs.A
#finding the first eigenvalue (finite rate of increase)
dom.pos <- which.max(eigs.A[["values"]])
L1mean <- Re(eigs.A[["values"]][dom.pos])
L1mean
#sensitivity of projection matrices 
vw.s <- v %*% t (w) 
S <- vw.s/as.numeric(v %*% w)
#elasticity of projection matrices 
elas <- (A/L1mean) * S 
elasticity <- round(elas, 3)
### figure 3 - plot the proportional sensitivity to changes in F, P and G 
stage <- c("E", "J", "A")
F <- c(elasticity[1, 1:3])
P <- c(elasticity[1,1], elasticity[2,2], elasticity[3,3])
G <- c(elasticity[2,1], elasticity[3,2], NA)
sensitivities <- data.frame(stage, F, P, G)
sensitivity <- read.csv("cheat for now.csv") 
sens <- gather(sensitivities, vr, elas, F, P, G)
#change sensitivites data frame into the correct format 
fig.3 <- ggplot(sens, aes(x = stage, y = elas, colour = vr, vr)) + geom_line() + geom_point(size = 4) + labs(x = "Stage", y = "Elasticity")
fig.3 + scale_x_discrete(limits=c("E","J","A")) 

##### figure 2 
# Looking at the effect of Age of First Reproduction on Intrinsic rate of Increase (r)
#-------------- Calculating changes in rate of increase r resulting from subtracting and adding 1 yr to the 
#calculations of P, and G, for each of the three immature stages.
yellow_di1 <- mutate(yellow, di = ifelse(stage == "2", di + 1, di * 1))
yellow_dii1 <- mutate(yellow_di1, di = ifelse(stage == "3", di - 1, di * 1)) # stage 2 as increasing age at which first breeding is = increasing number of years as juvenile
yellow_di2 <- mutate(yellow, di = ifelse(stage == "2", di + 2, di * 1)) # ^therefore juvenile = birds which have fledged the nest but have not began attempts at breeding
yellow_dii2 <- mutate(yellow_di2, di = ifelse(stage == "3", di - 2, di * 1)) 
yellow_di3 <- mutate(yellow, di = ifelse(stage == "2", di + 3, di * 1)) #dii1 = di is stage duration, i1 means increased by 1
yellow_dii3 <- mutate(yellow_di3, di = ifelse(stage == "3", di - 3, di * 1))


#applying my matrix function to the 3 new tables 
matinc1 <- ysameanFunc(yellow_dii1)
matinc2 <- ysameanFunc(yellow_dii2)
matinc3 <- ysameanFunc(yellow_dii3)

#eigenanalyses 
eigs.matinc1 <- eigen(matinc1)
eigs.matinc2 <- eigen(matinc2)
eigs.matinc3 <- eigen(matinc3)


#finding the first eigenvalue (finite rate of increase)
dom.pos.i1 <- which.max(eigs.matinc1[["values"]])
L1 <- Re(eigs.matinc1[["values"]][dom.pos.i1])
dom.pos.i2 <- which.max(eigs.matinc2[["values"]])
L2 <- Re(eigs.matinc2[["values"]][dom.pos.i2])
dom.pos.i3 <- which.max(eigs.matinc3[["values"]])
L3 <- Re(eigs.matinc3[["values"]][dom.pos.i3])

#finding r
r1 <- log(L1)
r2 <- log(L2)
r3 <- log(L3)


#plotting 
age <- c(2, 3, 4, 5)
lambdas<- c(L1mean, L1, L2, L3)
rs <- log(lambdas)
table_agerep <- data.frame(age, rs)
graph <- ggplot(table_agerep, aes(x = age, y = rs)) + geom_line(size=1) + geom_point(size=2)
graph1 <- graph + labs(x = "Age of First Reproduction (yr)", y = "Intrinsic rate of Increase (r)")  
figure2 <- graph1 + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())
figure2

