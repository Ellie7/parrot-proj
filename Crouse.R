#Crouse et al. 1987 
#A STAGE-BASED POPULATION MODEL FOR LOGGERHEAD SEA TURTLES AND IMPLICATIONS FOR CONSERVATION
library(dplyr)
library(ggplot2) 
library(knitr)
library(primer)
table.3 <- read.csv("~/1 UNIVERSITY/Level 4/Project & Dissertation/Crouse 1987/table 3 from crouse.csv")
View(table.3)
parrot-proj/table 3 from crouse.csv
table.3
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
### making a function which creates the matrix 
lifetable <- table.3
myFunc <- function (lifetable) 
{
  fecs <- select(lifetable, fecundity)
  pi <- select(lifetable, annual_survivorship)
  di <-select(lifetable, stage_duration)
  Pi <- ((1 - (pi^(di - 1)))/(1 - (pi^di)))*pi
  Gi <- (pi^di*(1 - pi))/(1 - pi^di)
  mat1 <- matrix(0, nrow = 7, ncol = 7)
  for (i in 2:7) {
    for (j in 2:7) mat1[i, j]
  }
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

#recreating table 4
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
ggplot()
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
#sensitivity of projection matrices 
vw.s <- v %*% t (w) 
S <- (S <- vw.s/as.numeric(v %*% w)) 
#elasticity of projection matrices 
elas <- (A/L1) * S 
elasticity <- round(elas, 3)
### figure 3 - plot the proportional sensitivity to changes in F, P and G 
stage <- c(1:7)
F <- c(elasticity[1, 1:7])
P <- c(elasticity[1,1], elasticity[2,2], elasticity[3,3], elasticity[4,4], elasticity[5,5], elasticity[6,6], elasticity[7,7])
G <- c(elasticity[2,1], elasticity[3,2], elasticity[4,3], elasticity[5,4], elasticity[6,5], elasticity[7,6], 0)
sensitvities <- data.frame(stage, F, P, G)
sens <- read.csv("~/1 UNIVERSITY/Level 4/Project & Dissertation/Crouse 1987/sens.csv")
fig.3 <- ggplot(sens, aes(x = stage, y = sens, colour = supp, shape = supp)) + geom_line() + geom_point(size = 4) + labs(x = "Stage", y = "Elasticity")
fig.3 
#-------------- Calculating changes in rate of increase r resulting from simulated changes in fecundity and survival of individual life history stages in the loggerhead population matrix 
### figure 1 
#changes in rate of increase r resulting from simulated changes in fecundity and survival of individual life history 
#caculating r determined in the baseline run of the matrix
table.3.50 <- read.csv("~/1 UNIVERSITY/Level 4/Project & Dissertation/Crouse 1987/table 3 surv fec altered .csv")
View(table.3.50)
#first recalculate eigenvalues for a 50% decrease in survivorship & 50% decrease in fecundity
lifetable <- table.3.50
A <- myFunc(lifetable)
eigs.A <- eigen(A)
eigs.A
#finding the first eigenvalue (finite rate of increase)
dom.pos <- which.max(eigs.A[["values"]])
L1N <- Re(eigs.A[["values"]][dom.pos]) #N for survivorship to distibuish from initial matrix 
L1N #=0.4046335  
lambda <- eigs.A[["values"]]
exp <- (Re(eigs.A[["values"]]))^2
rs <- log(sqrt(exp))




stage <- c("Eggs/Hatchlings", "Small Juveniles", "Large Juveniles", "Subadults", "Novice Breeders", "1st-yr Remigrants", "Mature Breeders")
changes <- data.frame(stage, rs)
ggplot(changes, aes(x = stage, y = rs)) + geom_bar()