#Crouse et al. 1987 
#A STAGE-BASED POPULATION MODEL FOR LOGGERHEAD SEA TURTLES AND IMPLICATIONS FOR CONSERVATION
library(dplyr)
library(ggplot2) 
library(knitr)
library(primer)
library(Rmisc)
library(agricolae)
table.3 <- read.csv("table 3 from crouse.csv")
#----------------- Creating a stage-based projection matrix, for each stage, 
#calculating the repro- ductive output (F,), the probability of surviving and 
#growing into the next stage (G,). and the probability of surviving and 
#remaining in the same stage (P,).
#Fs
fecs <- select(table.3, fecundity)
#starting values 
pi <- select(table.3, annual_survivorship)
di <-select(table.3, stage_duration)
#Ps
Pi <- ((1 - (pi^(di - 1)))/(1 - (pi^di)))*pi
# Pi's
#1              0.0000
#2              0.7371
#3              0.6611
#4              0.6907
#5              0.0000
#6              0.0000
#7              0.8086

#Gs
Gi <- (pi^di*(1 - pi))/(1 - pi^di) 
#life table (parts of anyway)
life_table <- data.frame(fecs, Gi, Pi)
##making the matrix: 
mat1 <- matrix(NA, nrow = 7, ncol = 7, byrow = T) 
#fill with fecundity 
mat1[1,] <- life_table$fecundity 
#add Ps 
mat1[1,1] <- Pi$annual_survivorship[1] 
mat1[2,2] <- Pi$annual_survivorship[2] 
mat1[3,3] <- Pi$annual_survivorship[3]
mat1[4,4] <- Pi$annual_survivorship[4]
mat1[5,5] <- Pi$annual_survivorship[5]
mat1[6,6] <- Pi$annual_survivorship[6]
mat1[7,7] <- Pi$annual_survivorship[7]
mat1
#add Gs 
mat1[2,1] <- Gi$annual_survivorship[1]
mat1[3,2] <- Gi$annual_survivorship[2]
mat1[4,3] <- Gi$annual_survivorship[3]
mat1[5,4] <- Gi$annual_survivorship[4]
mat1[6,5] <- Gi$annual_survivorship[5]
mat1[7,6] <- Gi$annual_survivorship[6] 
mat1
#remove NAs
#shows location of NAs
is.na(mat1)
#replaces NAs with 0s
mat1[is.na(mat1)] <- 0
mat1#almost there 
A <- mat1
### making a function which creates the matrix (like the code above just in the form of a function)
lifetable <- table.3
myFunc <- function (lifetable) 
{
  fecs <- select(lifetable, fecundity)
  pi <- select(lifetable, annual_survivorship)
  di <-select(lifetable, stage_duration)
  Pi <- ((1 - (pi^(di - 1)))/(1 - (pi^di)))*pi
  Gi <- (pi^di*(1 - pi))/(1 - pi^di)
  mat1 <- matrix(0, nrow = 7, ncol = 7)
    #add Fs
  mat1[1,] <- lifetable$fecundity 
  #add Ps (diagonals)
  mat1[1,1] <- Pi$annual_survivorship[1] 
  mat1[2,2] <- Pi$annual_survivorship[2] 
  mat1[3,3] <- Pi$annual_survivorship[3]
  mat1[4,4] <- Pi$annual_survivorship[4]
  mat1[5,5] <- Pi$annual_survivorship[5]
  mat1[6,6] <- Pi$annual_survivorship[6]
  mat1[7,7] <- Pi$annual_survivorship[7]
  mat1
  #add Gs (off-diagonals)
  mat1[2,1] <- Gi$annual_survivorship[1]
  mat1[3,2] <- Gi$annual_survivorship[2]
  mat1[4,3] <- Gi$annual_survivorship[3]
  mat1[5,4] <- Gi$annual_survivorship[4]
  mat1[6,5] <- Gi$annual_survivorship[5]
  mat1[7,6] <- Gi$annual_survivorship[6] 
  return(mat1)
}

myFunc(lifetable) 

# manually recreating table 4 (produces the same matrix as the code above just manually)
B<- matrix(c(0, 0, 0, 0, 127, 4, 80, 0.6747, 0.7370, 0, 0,0, 0, 0, 0, 0.0486, 0.6610, 0, 0, 
              0, 0, 0, 0, 0.0147, 0.6907, 0, 0, 0, 0, 0, 0, 0.0518, 0, 0, 0, 0, 0, 0, 0, 0.8091, 
              0, 0, 0, 0, 0, 0, 0, 0.8091, 0.8089), nr=7, byrow = TRUE) 

#--------------- population projections 
#stage structure growth (multiple steps)
N0 <- matrix(c(10000,10000,10000,10000,10000,10000,10000), ncol=1)
years <- 20
N.projections <- matrix(0, nrow = nrow(A), ncol = years + 1) 
N.projections[,1] <- N0 
for (i in 1:years) 
  {N.projections[, i + 1] <- A %*% N.projections[, i] 
matplot(0:years, t(N.projections), type = "l", lty = 1:3, 
        col = 1, ylab = "Stage Abundance", xlab = "Year")}
#annual growth rate
N.totals <- apply(N.projections, 2, sum)
Rs <- N.totals[-1]/N.totals[-(years + 1)]
#eigen analysis 
eigs.A <- eigen(A)
eigs.A
#finding the first eigenvalue (finite rate of increase)
dom.pos <- which.max(eigs.A[["values"]])
L1 <- Re(eigs.A[["values"]][dom.pos])
L1
#=0.9451619
#finding r 
r <- log(L1)
r
# r = -0.056399
#power method
t <- 20
Nt <- N0/sum(N0)
R.t <- numeric(t)
for (i in 1:t) R.t[i] <- {
  Nt1 <- A %*% Nt
  R <- sum(Nt1)/sum(Nt)
  R
} 
#You might need to adjust the number of iterations to make sure the
#value has stabilised (how can you tell that it has?).

#calculating the stable stage distribution 
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
#1.000000   1.400863   5.997875 115.559274 567.386379 505.836132 585.956078
#to create table 5 (from Crouse 1987)
tab <- select(table.3, stage_number, class)
tab_5<- data.frame(tab, stable, RV)
colnames(tab_5) <- c("Stage number", "Stage Class", "Stable stage distribution (Dominant eigenvector)", "Reproductive values (left eigenvector)" )
tab_5
kable(tab_5, caption = "Table 5. Stable stage distribution (wJ) and reproductive values (v') for the loggerhead population matrix given in Table 4.")

#---------------------------------- sensitivity analyses 

#figure 1
#-------------- Calculating changes in rate of increase r resulting from simulated changes in fecundity and survival of individual life history 
#stages in the loggerhead population matrix 
#decresing fecundity and survival by 50%
table.3_fecAdjust<-mutate(table.3, fecundity = 0.5*fecundity)
table.3_Surv1Adj <- mutate(table.3, annual_survivorship = ifelse(stage_number == "1", annual_survivorship * 0.5, annual_survivorship * 1))
table.3_Surv2Adj <- mutate(table.3, annual_survivorship = ifelse(stage_number == "2", annual_survivorship * 0.5, annual_survivorship * 1))
table.3_Surv3Adj <- mutate(table.3, annual_survivorship = ifelse(stage_number == "3", annual_survivorship * 0.5, annual_survivorship * 1))
table.3_Surv4Adj <- mutate(table.3, annual_survivorship = ifelse(stage_number == "4", annual_survivorship * 0.5, annual_survivorship * 1))
table.3_Surv5Adj <- mutate(table.3, annual_survivorship = ifelse(stage_number == "5", annual_survivorship * 0.5, annual_survivorship * 1))
table.3_Surv6Adj <- mutate(table.3, annual_survivorship = ifelse(stage_number == "6", annual_survivorship * 0.5, annual_survivorship * 1))
table.3_Surv7Adj <- mutate(table.3, annual_survivorship = ifelse(stage_number == "7", annual_survivorship * 0.5, annual_survivorship * 1))
#50% increase in fecundity or an increase in survivorship to 1.0.
table.3_fecInc<-mutate(table.3, fecundity = 1.5*fecundity)
table.3_Surv1Inc <- mutate(table.3, annual_survivorship = ifelse(stage_number == "1", (0.999-annual_survivorship+annual_survivorship), annual_survivorship * 1)) # 0.999 to test if it's a problem with 1
table.3_Surv2Inc <- mutate(table.3, annual_survivorship = ifelse(stage_number == "2", (0.999-annual_survivorship+annual_survivorship), annual_survivorship * 1))
table.3_Surv3Inc <- mutate(table.3, annual_survivorship = ifelse(stage_number == "3", (0.999-annual_survivorship+annual_survivorship), annual_survivorship * 1))
table.3_Surv4Inc <- mutate(table.3, annual_survivorship = ifelse(stage_number == "4", (0.999-annual_survivorship+annual_survivorship), annual_survivorship * 1))
table.3_Surv5Inc <- mutate(table.3, annual_survivorship = ifelse(stage_number == "5", (0.999-annual_survivorship+annual_survivorship), annual_survivorship * 1))
table.3_Surv6Inc <- mutate(table.3, annual_survivorship = ifelse(stage_number == "6", (0.999-annual_survivorship+annual_survivorship), annual_survivorship * 1))
table.3_Surv7Inc <- mutate(table.3, annual_survivorship = ifelse(stage_number == "7", (0.999-annual_survivorship+annual_survivorship), annual_survivorship * 1))
#eigenanalysis 
#decreases 
AFd <- myFunc(table.3_fecAdjust) 
eigs.AFd <- eigen(AFd)
A1d <- myFunc(table.3_Surv1Adj) 
eigs.A1d <- eigen(A1d)
A2d <- myFunc(table.3_Surv2Adj) 
eigs.A2d <- eigen(A2d)
A3d <- myFunc(table.3_Surv3Adj) 
eigs.A3d <- eigen(A3d)
A4d <- myFunc(table.3_Surv4Adj) 
eigs.A4d <- eigen(A4d)
A5d <- myFunc(table.3_Surv5Adj) 
eigs.A5d <- eigen(A5d)
A6d <- myFunc(table.3_Surv6Adj) 
eigs.A6d <- eigen(A6d)
A7d <- myFunc(table.3_Surv7Adj) 
eigs.A7d <- eigen(A7d)
#increases 
AFi <- myFunc(table.3_fecInc) # F for fecundity, i for increase (50%)
eigs.AFi <- eigen(AFi)
A1i <- myFunc(table.3_Surv1Inc) # 1 for stage number, i for increase (to 1.0)
eigs.A1i <- eigen(A1i)
A2i <- myFunc(table.3_Surv2Inc) 
eigs.A2i <- eigen(A2i)
A3i <- myFunc(table.3_Surv3Inc) 
eigs.A3i <- eigen(A3i)
A4i <- myFunc(table.3_Surv4Inc) 
eigs.A4i <- eigen(A4i)
A5i <- myFunc(table.3_Surv5Inc) 
eigs.A5i <- eigen(A5i)
A6i <- myFunc(table.3_Surv6Inc) 
eigs.A6i <- eigen(A6i)
A7i <- myFunc(table.3_Surv7Inc) 
eigs.A7i <- eigen(A7i)
#finding the first eigenvalue (finite rate of increase)
#decreases
dom.pos.Fd <- which.max(eigs.AFd[["values"]])
LFd <- Re(eigs.AFd[["values"]][dom.pos.Fd]) 
dom.pos.1d <- which.max(eigs.A1d[["values"]])
L1d <- Re(eigs.A1d[["values"]][dom.pos.1d]) 
dom.pos.2d <- which.max(eigs.A2d[["values"]])
L2d <- Re(eigs.A2d[["values"]][dom.pos.2d]) 
dom.pos.3d <- which.max(eigs.A3d[["values"]])
L3d <- Re(eigs.A3d[["values"]][dom.pos.3d]) 
dom.pos.4d <- which.max(eigs.A4d[["values"]])
L4d <- Re(eigs.A4d[["values"]][dom.pos.4d]) 
dom.pos.5d <- which.max(eigs.A5d[["values"]])
L5d <- Re(eigs.A5d[["values"]][dom.pos.5d]) 
dom.pos.6d <- which.max(eigs.A6d[["values"]])
L6d <- Re(eigs.A6d[["values"]][dom.pos.6d]) 
dom.pos.7d <- which.max(eigs.A7d[["values"]])
L7d <- Re(eigs.A7d[["values"]][dom.pos.7d]) 
#increases 
dom.pos.Fi <- which.max(eigs.AFi[["values"]])
LFi <- Re(eigs.AFi[["values"]][dom.pos.Fi]) 
dom.pos.1i <- which.max(eigs.A1i[["values"]])
L1i <- Re(eigs.A1i[["values"]][dom.pos.1i])
dom.pos.2i <- which.max(eigs.A2i[["values"]])
L2i <- Re(eigs.A2i[["values"]][dom.pos.2i]) 
dom.pos.3i <- which.max(eigs.A3i[["values"]])
L3i <- Re(eigs.A3i[["values"]][dom.pos.3i]) 
dom.pos.4i <- which.max(eigs.A4i[["values"]])
L4i <- Re(eigs.A4i[["values"]][dom.pos.4i]) 
dom.pos.5i <- which.max(eigs.A5i[["values"]])
L5i <- Re(eigs.A5i[["values"]][dom.pos.5i]) 
dom.pos.6i <- which.max(eigs.A6i[["values"]])
L6i <- Re(eigs.A6i[["values"]][dom.pos.6i]) 
dom.pos.7i <- which.max(eigs.A7i[["values"]])
L7i <- Re(eigs.A7i[["values"]][dom.pos.7i]) 

#plotting figure 1
stage_class <- c("Fecundity", "Eggs/Hatchlings", "Small Juveniles", "Large Juveniles", "Subadults", "Novice Breeders", "1st-yr Remigrants", "Mature Breeders")
lambdas<- c(LFd, L1d, L2d, L3d, L4d, L5d, L6d, L7d)
rs <- log(lambdas)
table_decrease <- data.frame(stage_class, rs)
graph <- ggplot(table_decrease, aes(x = stage_class, y = rs)) + geom_bar(stat = "identity")
graph1 <- graph + labs(x = "Stage Class", y = "Intrinsic rate of Increase (r)") + geom_hline(aes(yintercept=-0.056399), linetype = "dashed", size = 0.75)  
graph2 <- graph1 + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ expand_limits(y = -0.25, x = 0.2)
figure1a <- graph2 + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + scale_x_discrete(limits = c("Fecundity", "Eggs/Hatchlings", "Small Juveniles", "Large Juveniles", "Subadults", "Novice Breeders", "1st-yr Remigrants", "Mature Breeders"))
figure1a 
stage_class <- c("Fecundity", "Eggs/Hatchlings", "Small Juveniles", "Large Juveniles", "Subadults", "Novice Breeders", "1st-yr Remigrants", "Mature Breeders")
lambdas<- c(LFi, L1i, L2i, L3i, L4i, L5i, L6i, L7i)
rs <- log(lambdas)
table_decrease <- data.frame(stage_class, rs)
graph <- ggplot(table_decrease, aes(x = stage_class, y = rs)) + geom_bar(stat = "identity")
graph1 <- graph + labs(x = "Stage Class", y = "Intrinsic rate of Increase (r)") + geom_hline(aes(yintercept=-0.056399), linetype = "dashed", size = 0.75)  
graph2 <- graph1 + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ expand_limits(y = -0.25, x = 0.2)
figure1b <- graph2 + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + scale_x_discrete(limits = c("Fecundity", "Eggs/Hatchlings", "Small Juveniles", "Large Juveniles", "Subadults", "Novice Breeders", "1st-yr Remigrants", "Mature Breeders"))
figure1b 

##### figure 2 
#These conditions, a 6-yr decrease or increase in tthe age of reproductive maturity, were simulated (Fig. 2) by subtracting and adding 2 yr to the 
#calculations of P, and G, for each of the three immature stages. 
# In fact, a mere 3-yr reduction in the age of first repro- duction, well within the bounds of the growth estimates available, 
#comes very close to halting the decline in this population
#-------------- Calculating changes in rate of increase r resulting from subtracting and adding 2 yr to the 
#calculations of P, and G, for each of the three immature stages.
table.3_i6 <- mutate(table.3, stage_duration = ifelse(stage_number == "2", stage_duration + 2, stage_duration * 1))
table.3_in6 <- mutate(table.3_i6, stage_duration = ifelse(stage_number == "3", stage_duration + 2, stage_duration * 1))#for some reasion had to do each separately for it work
table.3_inc6 <- mutate(table.3_in6,stage_duration = ifelse(stage_number =="4", stage_duration + 2, stage_duration * 1))#final table for increased by 6 
table.3_d6 <- mutate(table.3, stage_duration = ifelse(stage_number == "2", stage_duration - 2, stage_duration * 1))
table.3_de6 <- mutate(table.3_d6, stage_duration = ifelse(stage_number == "3", stage_duration - 2, stage_duration * 1))
table.3_dec6 <- mutate(table.3_de6,stage_duration = ifelse(stage_number =="4", stage_duration - 2, stage_duration * 1))
table.3_d3 <- mutate(table.3, stage_duration = ifelse(stage_number == "2", stage_duration - 1, stage_duration * 1))
table.3_de3 <- mutate(table.3_d3, stage_duration = ifelse(stage_number == "3", stage_duration - 1, stage_duration * 1))
table.3_dec3 <- mutate(table.3_de3,stage_duration = ifelse(stage_number =="4", stage_duration - 1, stage_duration * 1))
#applying my matrix function to the 3 new tables 
mat28 <- myFunc(table.3_inc6) 
mat16 <- myFunc(table.3_dec6) 
mat19 <- myFunc(table.3_dec3)
#eigenanalyses 
eigs.mat28 <- eigen(mat28)
eigs.mat16 <- eigen(mat16)
eigs.mat19 <- eigen(mat19)
#finding the first eigenvalue (finite rate of increase)
dom.pos28 <- which.max(eigs.mat28[["values"]])
L28 <- Re(eigs.mat28[["values"]][dom.pos28])
dom.pos16 <- which.max(eigs.mat16[["values"]])
L16 <- Re(eigs.mat16[["values"]][dom.pos16])
dom.pos19 <- which.max(eigs.mat19[["values"]])
L19 <- Re(eigs.mat19[["values"]][dom.pos19])
#finding r
r28 <- log(L28)
r16 <- log(L16)
r19 <- log(L19)
#plotting 
age <- c(16, 19, 22, 28)
lambdas<- c(L16, L19, L1, L28)
rs <- log(lambdas)
table_agerep <- data.frame(age, rs)
graph <- ggplot(table_agerep, aes(x = age, y = rs)) + geom_line(size=1) + geom_point(size=2)
graph1 <- graph + labs(x = "Age of First Reproduction (yr)", y = "Intrinsic rate of Increase (r)") + geom_hline(aes(yintercept=0), size = 0.5)  
figure2 <- graph1 + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())
figure2

######## figure 3 
#the elasticity, or proportional sensitivity of lambda to changes in fecundity F, survival while remaining in the same stage P, and survival 
# with growth, G. Because the elasticities of these matrix elements sum to 1, they can be compared directly in terms of their contribution to the 
# population growth rate 
#sensitivity of projection matrices 
vw.s <- v %*% t (w) 
S <- vw.s/as.numeric(v %*% w)
#elasticity of projection matrices 
elas <- (A/L1) * S 
elasticity <- round(elas, 3)
### figure 3 - plot the proportional sensitivity to changes in F, P and G 
stage <- c(1:7)
F <- c(elasticity[1, 1:7])
P <- c(elasticity[1,1], elasticity[2,2], elasticity[3,3], elasticity[4,4], elasticity[5,5], elasticity[6,6], elasticity[7,7])
G <- c(elasticity[2,1], elasticity[3,2], elasticity[4,3], elasticity[5,4], elasticity[6,5], elasticity[7,6], 0)
sensitvities <- data.frame(stage, F, P, G)
sens <- read.csv("~/1 UNIVERSITY/Level 4/Project & Dissertation/Crouse 1987/sens.csv") #been lazy and just used excel rather than finding how to 
#change sensitivites data frame into the correct format 
fig.3 <- ggplot(sens, aes(x = stage, y = sens, colour = supp, shape = supp)) + geom_line() + geom_point(size = 4) + labs(x = "Stage", y = "Elasticity")
fig.3 

############# figure 4 
#(a) The elasticity, or proportional sensitivity, of XA, to changes in annual stage-specific survival probability pi. (b) The elasticity of lambda to 
#changes in stage duration di. Elasticity in stage duration is negative because stage duration and population growth rates r are inversely related.
#sensitivity of projection matrices 
#In this model, P, and G, are derived parameters; they depend on both stage-specific annual survival probability p and stage duration d. 
#We have also calculated the elasticities of lambda with respect to these parameters
#change p values
life_pi1 <- mutate(lifetable, annual_survivorship = ifelse(stage_number == "1", annual_survivorship * 1.1, annual_survivorship * 1)) # p for survival probility, i for increase, 1 for the stage that its p is being altered
life_pd1 <- mutate(lifetable, annual_survivorship = ifelse(stage_number == "1", annual_survivorship * 0.9, annual_survivorship * 1)) # d for decrease
life_pi2 <- mutate(lifetable, annual_survivorship = ifelse(stage_number == "2", annual_survivorship * 1.1, annual_survivorship * 1))
life_pd2 <- mutate(lifetable, annual_survivorship = ifelse(stage_number == "2", annual_survivorship * 0.9, annual_survivorship * 1))
life_pi3 <- mutate(lifetable, annual_survivorship = ifelse(stage_number == "3", annual_survivorship * 1.1, annual_survivorship * 1))
life_pd3 <- mutate(lifetable, annual_survivorship = ifelse(stage_number == "3", annual_survivorship * 0.9, annual_survivorship * 1))
life_pi4 <- mutate(lifetable, annual_survivorship = ifelse(stage_number == "4", annual_survivorship * 1.1, annual_survivorship * 1))
life_pd4 <- mutate(lifetable, annual_survivorship = ifelse(stage_number == "4", annual_survivorship * 0.9, annual_survivorship * 1))
life_pi5 <- mutate(lifetable, annual_survivorship = ifelse(stage_number == "5", annual_survivorship * 1.1, annual_survivorship * 1))
life_pd5 <- mutate(lifetable, annual_survivorship = ifelse(stage_number == "5", annual_survivorship * 0.9, annual_survivorship * 1))
life_pi6 <- mutate(lifetable, annual_survivorship = ifelse(stage_number == "6", annual_survivorship * 1.1, annual_survivorship * 1))
life_pd6 <- mutate(lifetable, annual_survivorship = ifelse(stage_number == "6", annual_survivorship * 0.9, annual_survivorship * 1))
life_pi7 <- mutate(lifetable, annual_survivorship = ifelse(stage_number == "7", annual_survivorship * 1.1, annual_survivorship * 1))
life_pd7 <- mutate(lifetable, annual_survivorship = ifelse(stage_number == "7", annual_survivorship * 0.9, annual_survivorship * 1))
#ds increasing and decreasing by 5% 
life_di1 <- mutate(lifetable, stage_duration = ifelse(stage_number == "1", annual_survivorship * 1.1, annual_survivorship * 1))
life_dd1 <- mutate(lifetable, stage_duration = ifelse(stage_number == "1", annual_survivorship * 0.9, annual_survivorship * 1))
life_di2 <- mutate(lifetable, stage_duration = ifelse(stage_number == "2", annual_survivorship * 1.1, annual_survivorship * 1))
life_dd2 <- mutate(lifetable, stage_duration = ifelse(stage_number == "2", annual_survivorship * 0.9, annual_survivorship * 1))
life_di3 <- mutate(lifetable, stage_duration = ifelse(stage_number == "3", annual_survivorship * 1.1, annual_survivorship * 1))
life_dd3 <- mutate(lifetable, stage_duration = ifelse(stage_number == "3", annual_survivorship * 0.9, annual_survivorship * 1))
life_di4 <- mutate(lifetable, stage_duration = ifelse(stage_number == "4", annual_survivorship * 1.1, annual_survivorship * 1))
life_dd4 <- mutate(lifetable, stage_duration = ifelse(stage_number == "4", annual_survivorship * 0.9, annual_survivorship * 1))
life_di5 <- mutate(lifetable, stage_duration = ifelse(stage_number == "5", annual_survivorship * 1.1, annual_survivorship * 1))
life_dd5 <- mutate(lifetable, stage_duration = ifelse(stage_number == "5", annual_survivorship * 0.9, annual_survivorship * 1))
life_di6 <- mutate(lifetable, stage_duration = ifelse(stage_number == "6", annual_survivorship * 1.1, annual_survivorship * 1))
life_dd6 <- mutate(lifetable, stage_duration = ifelse(stage_number == "6", annual_survivorship * 0.9, annual_survivorship * 1))
life_di7 <- mutate(lifetable, stage_duration = ifelse(stage_number == "7", annual_survivorship * 1.1, annual_survivorship * 1))
life_dd7 <- mutate(lifetable, stage_duration = ifelse(stage_number == "7", annual_survivorship * 0.9, annual_survivorship * 1))
#make new matrices
mat_pi1 <- myFunc(life_pi1)
mat_pd1 <- myFunc(life_pd1)
mat_pi2 <- myFunc(life_pi2)
mat_pd2 <- myFunc(life_pd2)
mat_pi3 <- myFunc(life_pi3)
mat_pd3 <- myFunc(life_pd3)
mat_pi4 <- myFunc(life_pi4)
mat_pd4 <- myFunc(life_pd4)
mat_pi5 <- myFunc(life_pi5)
mat_pd5 <- myFunc(life_pd5)
mat_pi6<- myFunc(life_pi6)
mat_pd6 <- myFunc(life_pd6)
mat_pi7 <- myFunc(life_pi7)
mat_pd7 <- myFunc(life_pd7) 
mat_di1 <- myFunc(life_di1)
mat_dd1 <- myFunc(life_dd1)
mat_di2 <- myFunc(life_di2)
mat_dd2 <- myFunc(life_dd2)
mat_di3 <- myFunc(life_di3)
mat_dd3 <- myFunc(life_dd3)
mat_di4 <- myFunc(life_di4)
mat_dd4 <- myFunc(life_dd4)
mat_di5 <- myFunc(life_di5)
mat_dd5 <- myFunc(life_dd5)
mat_di6 <- myFunc(life_di6)
mat_dd6 <- myFunc(life_dd6)
mat_di7 <- myFunc(life_di7)
mat_dd7 <- myFunc(life_dd7)
#sensitivity of projection matrices
#eigen analysis 
eigs_pi1 <- eigen(mat_pi1)
eigs_di1 <- eigen(mat_di1)
eigs_pi2 <- eigen(mat_pi2)
eigs_di2 <- eigen(mat_di2)
eigs_pi3 <- eigen(mat_pi3)
eigs_di3 <- eigen(mat_di3)
eigs_pi4 <- eigen(mat_pi4)
eigs_di4 <- eigen(mat_di4)
eigs_pi5 <- eigen(mat_pi5)
eigs_di5 <- eigen(mat_di5)
eigs_pi6 <- eigen(mat_pi6)
eigs_di6 <- eigen(mat_di6)
eigs_pi7 <- eigen(mat_pi7)
eigs_di7 <- eigen(mat_di7)
eigs_pd1 <- eigen(mat_pd1)
eigs_dd1 <- eigen(mat_dd1)
eigs_pd2 <- eigen(mat_pd2)
eigs_dd2 <- eigen(mat_dd2)
eigs_pd3 <- eigen(mat_pd3)
eigs_dd3 <- eigen(mat_dd3)
eigs_pd4 <- eigen(mat_pd4)
eigs_dd4 <- eigen(mat_dd4)
eigs_pd5 <- eigen(mat_pd5)
eigs_dd5 <- eigen(mat_dd5)
eigs_pd6 <- eigen(mat_pd6)
eigs_dd6 <- eigen(mat_dd6)
eigs_pd7 <- eigen(mat_pd7)
eigs_dd7 <- eigen(mat_dd7)
#finding first eigenvalue
#1
dom.pos.pi1 <- which.max(eigs_pi1[["values"]])
Lpi1 <- Re(eigs_pi1[["values"]][dom.pos.pi1]) 
dom.pos.di1 <- which.max(eigs_di1[["values"]])
Ldi1 <- Re(eigs_di1[["values"]][dom.pos.di1])
dom.pos.pd1 <- which.max(eigs_pd1[["values"]])
Lpd1 <- Re(eigs_pd1[["values"]][dom.pos.pd1]) 
dom.pos.dd1 <- which.max(eigs_dd1[["values"]])
Ldd1 <- Re(eigs_dd1[["values"]][dom.pos.dd1])
#2
dom.pos.pi2 <- which.max(eigs_pi2[["values"]])
Lpi2 <- Re(eigs_pi2[["values"]][dom.pos.pi2]) 
dom.pos.di2 <- which.max(eigs_di2[["values"]])
Ldi2 <- Re(eigs_di2[["values"]][dom.pos.di2])
dom.pos.pd2 <- which.max(eigs_pd2[["values"]])
Lpd2 <- Re(eigs_pd2[["values"]][dom.pos.pd2]) 
dom.pos.dd2 <- which.max(eigs_dd2[["values"]])
Ldd2 <- Re(eigs_dd2[["values"]][dom.pos.dd2])
#3
dom.pos.pi3 <- which.max(eigs_pi3[["values"]])
Lpi3 <- Re(eigs_pi3[["values"]][dom.pos.pi3]) 
dom.pos.di3 <- which.max(eigs_di3[["values"]])
Ldi3 <- Re(eigs_di3[["values"]][dom.pos.di3])
dom.pos.pd3 <- which.max(eigs_pd3[["values"]])
Lpd3 <- Re(eigs_pd3[["values"]][dom.pos.pd3]) 
dom.pos.dd3 <- which.max(eigs_dd3[["values"]])
Ldd3 <- Re(eigs_dd3[["values"]][dom.pos.dd3])
#4
dom.pos.pi4 <- which.max(eigs_pi4[["values"]])
Lpi4 <- Re(eigs_pi4[["values"]][dom.pos.pi4]) 
dom.pos.di4 <- which.max(eigs_di4[["values"]])
Ldi4 <- Re(eigs_di4[["values"]][dom.pos.di4])
dom.pos.pd4 <- which.max(eigs_pd4[["values"]])
Lpd4 <- Re(eigs_pd4[["values"]][dom.pos.pd4]) 
dom.pos.dd4 <- which.max(eigs_dd4[["values"]])
Ldd4 <- Re(eigs_dd4[["values"]][dom.pos.dd4])
#5
dom.pos.pi5 <- which.max(eigs_pi5[["values"]])
Lpi5 <- Re(eigs_pi5[["values"]][dom.pos.pi5]) 
dom.pos.di5 <- which.max(eigs_di5[["values"]])
Lpd5 <- Re(eigs_di5[["values"]][dom.pos.di5])
dom.pos.pd5 <- which.max(eigs_pd5[["values"]])
Lpd5 <- Re(eigs_pd5[["values"]][dom.pos.pd5]) 
dom.pos.dd5 <- which.max(eigs_dd5[["values"]])
Ldd5 <- Re(eigs_dd5[["values"]][dom.pos.dd5])
#6 
dom.pos.pi6 <- which.max(eigs_pi6[["values"]])
Lpi6 <- Re(eigs_pi6[["values"]][dom.pos.pi6]) 
dom.pos.di6 <- which.max(eigs_di6[["values"]])
Lpi6 <- Re(eigs_di6[["values"]][dom.pos.di6])
dom.pos.pd6 <- which.max(eigs_pd6[["values"]])
Lpd6 <- Re(eigs_pd6[["values"]][dom.pos.pd6]) 
dom.pos.dd6 <- which.max(eigs_dd6[["values"]])
Ldd6 <- Re(eigs_dd6[["values"]][dom.pos.dd6])
#7 
dom.pos.pi7 <- which.max(eigs_pi7[["values"]])
Lpi7 <- Re(eigs_pi7[["values"]][dom.pos.pi7]) 
dom.pos.di7 <- which.max(eigs_di7[["values"]])
Lpi7 <- Re(eigs_di7[["values"]][dom.pos.di7])
dom.pos.pd7 <- which.max(eigs_pd7[["values"]])
Lpd7 <- Re(eigs_pd7[["values"]][dom.pos.pd7]) 
dom.pos.dd7 <- which.max(eigs_dd7[["values"]])
Ldd7 <- Re(eigs_dd7[["values"]][dom.pos.dd7])

wd1 <- Re(eigs_pd1[["vectors"]][, dom.pos]) #just 5% decrese in p for now 
wd2 <- Re(eigs_pd2[["vectors"]][, dom.pos])
wd3 <- Re(eigs_pd3[["vectors"]][, dom.pos])
wd4 <- Re(eigs_pd4[["vectors"]][, dom.pos])
wd5 <- Re(eigs_pd5[["vectors"]][, dom.pos])
wd6 <- Re(eigs_pd6[["vectors"]][, dom.pos])
wd7 <- Re(eigs_pd7[["vectors"]][, dom.pos])

wdd1 <- Re(eigs_dd1[["vectors"]][, dom.pos]) #just 5% decrese in p for now 
wdd2 <- Re(eigs_dd2[["vectors"]][, dom.pos])
wdd3 <- Re(eigs_dd3[["vectors"]][, dom.pos])
wdd4 <- Re(eigs_dd4[["vectors"]][, dom.pos])
wdd5 <- Re(eigs_dd5[["vectors"]][, dom.pos])
wdd6 <- Re(eigs_dd6[["vectors"]][, dom.pos])
wdd7 <- Re(eigs_dd7[["vectors"]][, dom.pos])
#calculating reproductive value
t <- 20 #random for now 
#ps
Md1 <- eigen(t(mat_pd1))
Md2 <- eigen(t(mat_pd2))
Md3 <- eigen(t(mat_pd3))
Md4 <- eigen(t(mat_pd4))
Md5 <- eigen(t(mat_pd5))
Md6 <- eigen(t(mat_pd6))
Md7 <- eigen(t(mat_pd7))
#ds
Mdd1 <- eigen(t(mat_dd1))
Mdd2 <- eigen(t(mat_dd2))
Mdd3 <- eigen(t(mat_dd3))
Mdd4 <- eigen(t(mat_dd4))
Mdd5 <- eigen(t(mat_dd5))
Mdd6 <- eigen(t(mat_dd6))
Mdd7 <- eigen(t(mat_dd7))
#ps
vd1 <- Re(Md1$vectors[,which.max(Re(Md1$values))])
vd2 <- Re(Md2$vectors[,which.max(Re(Md2$values))])
vd3 <- Re(Md3$vectors[,which.max(Re(Md3$values))])
vd4 <- Re(Md4$vectors[,which.max(Re(Md4$values))])
vd5 <- Re(Md5$vectors[,which.max(Re(Md5$values))]) 
vd6 <- Re(Md6$vectors[,which.max(Re(Md6$values))])
vd7 <- Re(Md7$vectors[,which.max(Re(Md7$values))])
#ds 
vdd1 <- Re(Mdd1$vectors[,which.max(Re(Mdd1$values))])
vdd2 <- Re(Mdd2$vectors[,which.max(Re(Mdd2$values))])
vdd3 <- Re(Mdd3$vectors[,which.max(Re(Mdd3$values))])
vdd4 <- Re(Mdd4$vectors[,which.max(Re(Mdd4$values))])
vdd5 <- Re(Mdd5$vectors[,which.max(Re(Mdd5$values))]) 
vdd6 <- Re(Mdd6$vectors[,which.max(Re(Mdd6$values))])
vdd7 <- Re(Mdd7$vectors[,which.max(Re(Mdd7$values))])
#ps
vw.sd1 <- vd1 %*% t (wd1)
Sd1 <- vw.sd1/as.numeric(vd1 %*% wd1)
vw.sd2 <- vd2 %*% t (wd2)
Sd2 <- vw.sd2/as.numeric(vd2 %*% wd2)
vw.sd3 <- vd3 %*% t (wd3)
Sd3 <- vw.sd3/as.numeric(vd3 %*% wd3)
vw.sd4 <- vd4 %*% t (wd4)
Sd4 <- vw.sd4/as.numeric(vd4 %*% wd4)
vw.sd5 <- vd5 %*% t (wd5)
Sd5 <- vw.sd5/as.numeric(vd5 %*% wd5)
vw.sd6 <- vd6 %*% t (wd6)
Sd6 <- vw.sd6/as.numeric(vd6 %*% wd6)
vw.sd7 <- vd7 %*% t (wd7)
Sd7 <- vw.sd7/as.numeric(vd7 %*% wd7)
#ds 
vw.sdd1 <- vdd1 %*% t (wdd1)
Sdd1 <- vw.sdd1/as.numeric(vdd1 %*% wdd1)
vw.sdd2 <- vdd2 %*% t (wdd2)
Sdd2 <- vw.sdd2/as.numeric(vdd2 %*% wdd2)
vw.sdd3 <- vdd3 %*% t (wdd3)
Sdd3 <- vw.sdd3/as.numeric(vdd3 %*% wdd3)
vw.sdd4 <- vdd4 %*% t (wdd4)
Sdd4 <- vw.sd4/as.numeric(vdd4 %*% wdd4)
vw.sdd5 <- vdd5 %*% t (wdd5)
Sdd5 <- vw.sdd5/as.numeric(vdd5 %*% wdd5)
vw.sdd6 <- vdd6 %*% t (wdd6)
Sdd6 <- vw.sd6/as.numeric(vdd6 %*% wdd6)
vw.sdd7 <- vdd7 %*% t (wdd7)
Sdd7 <- vw.sdd7/as.numeric(vdd7 %*% wdd7)
#elasticity of projection matrices 
elasd1 <- (mat_pd1/Lpd1) * Sd1
elasticity1 <- round(elas, 3)
elasd2 <- (mat_pd2/Lpd1) * Sd2
elasticity2 <- round(elas, 3)
elasd3 <- (mat_pd3/Lpd1) * Sd3
elasticity3 <- round(elas, 3)
elasd4 <- (mat_pd4/Lpd1) * Sd4
elasticity4 <- round(elas, 3)
elasd5 <- (mat_pd5/Lpd1) * Sd5
elasticity5 <- round(elas, 3)
elasd6 <- (mat_pd6/Lpd1) * Sd6
elasticity6 <- round(elas, 3)
elasd7 <- (mat_pd7/Lpd1) * Sd7
elasticity7 <- round(elas, 3)

#elasticity of projection matrices for changes in d 
elasdd1 <- (mat_dd1/Ldd1) * Sdd1
elasticityd1 <- round(elas, 3)
elasdd2 <- (mat_dd2/Ldd1) * Sdd2
elasticityd2 <- round(elas, 3)
elasdd3 <- (mat_dd3/Ldd1) * Sdd3
elasticityd3 <- round(elas, 3)
elasdd4 <- (mat_dd4/Ldd1) * Sdd4
elasticityd4 <- round(elas, 3)
elasdd5 <- (mat_dd5/Ldd1) * Sdd5
elasticityd5 <- round(elas, 3)
elasdd6 <- (mat_dd6/Ldd1) * Sdd6
elasticityd6 <- round(elas, 3)
elasdd7 <- (mat_dd7/Ldd1) * Sdd7
elasticityd7 <- round(elas, 3)

### figure 4 - plot the proportional sensitivity to changes in survival probability, p, and stage duration, d. 
stage <- c(1:7)
P <- c(elasticity1[1,1], elasticity2[2,2], elasticity3[3,3], elasticity4[4,4], elasticity5[5,5], elasticity6[6,6], elasticity7[7,7])
sensi <- data.frame(stage, P)
fig.4a <- ggplot(sensi, aes(x = stage, y = P)) + geom_line() + geom_point(size = 4) + labs(x = "Stage", y = "Elasticity")
fig.4a 

stage <- c(1:7)
P <- c(elasticityd1[1,1], elasticityd2[2,2], elasticityd3[3,3], elasticityd4[4,4], elasticityd5[5,5], elasticityd6[6,6], elasticityd7[7,7])
sensi <- data.frame(stage, P)
fig.4b <- ggplot(sensi, aes(x = stage, y = P)) + geom_line() + geom_point(size = 4) + labs(x = "Stage", y = "Elasticity")
fig.4b 




