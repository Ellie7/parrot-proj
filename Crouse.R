#Crouse et al. 1987 
#A STAGE-BASED POPULATION MODEL FOR LOGGERHEAD SEA TURTLES AND IMPLICATIONS FOR CONSERVATION
install.packages("dplyr")
install.packages("ggplot2")
install.packages("knitr")
library(dplyr)
library(ggplot2) 
library(knitr)
table.3 <- read.csv("~/Parrot/table 3 from crouse.csv")
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
Pi <- round(Pi, 4)
#P1              0.0000
#P2              0.7198 (Crouse = 0.7370)
#P3              0.6535 (Crouse = 0.6610)
#P4              0.6675 (Crouse = 0.6907)
#P5              0.0000
#P6              0.0000
#P7              0.8086 (Crouse = 0.8089)

#Gs
Gi <- (pi^di*(1 - pi))/(1 - pi^di) 
Gi <- round(Gi, 4) 
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

#recreating table 4
A <- matrix(c(0, 0, 0, 0, 127, 4, 80, 0.6747, 0.7370, 0, 0,0, 0, 0, 0, 0.0486, 0.6610, 0, 0, 
              0, 0, 0, 0, 0.0147, 0.6907, 0, 0, 0, 0, 0, 0, 0.0518, 0, 0, 0, 0, 0, 0, 0, 0.8091, 
              0, 0, 0, 0, 0, 0, 0, 0.8091, 0.8089), nr=7, byrow = TRUE) 

#--------------- population projections 
#stage structure growth (multiple steps)
N0 <- matrix(c(1000,1000,1000,1000,1000,1000,1000), ncol=1)
years <- 10
N.projections <- matrix(0, nrow = nrow(A), ncol = years + 1) 
N.projections[,1] <- N0 
for (i in 1:years) 
  {N.projections[, i + 1] <- A %*% N.projections[, i] 
matplot(0:years, t(N.projections), type = "l", lty = 1:3, 
        col = 1, ylab = "Stage Abundance", xlab = "Year")}
#annual growth rate
N.totals <- apply(N.projections, 2, sum)
Rs <- N.totals[-1]/N.totals[-(years + 1)]
plot(0:(years - 1), Rs, type = "b", xlab = "Year", ylab = "R") 
#eigen analysis 
eigs.A <- eigen(A)
eigs.A
#finding the first eigenvalue (finite rate of increase)
dom.pos <- which.max(eigs.A[["values"]])
L1 <- Re(eigs.A[["values"]][dom.pos])
L1
#=0.945031
#power method
t <- 20
Nt <- N0/sum(N0)
R.t <- numeric(t)
for (i in 1:t) R.t[i] 
<- {
  Nt1 <- A %*% Nt
  R <- sum(Nt1)/sum(Nt)
  R
} 
#You might need to adjust the number of iterations to make sure the
#value has stabilised (how can you tell that it has?).
par(mar = c(5,4,3,2))
plot(1:t, R.t, type = "b", main = quote("Convergence Toward"* lambda))
points(t, L1, pch = 19, cex = 1.5)
#calculating the stable stage distribution 
w <- Re(eigs.A[["vectors"]][, dom.pos])
ssd <- w / sum(w)
stable <- round(ssd,3)
#0.207 0.670 0.115 0.007 0.000 0.000 0.002 
#calculating the reproductive value 
M <- eigen(t(A))
v <- Re(M$vectors[,which.max(Re(M$values))])
RV <- v / v[1]
RV 
#1.000000   1.400668   5.995523 115.844511 568.780852 507.373040 587.669314
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
plot(, Rs, type = "b", xlab = "Year", ylab = "R") 
#elasticity of projection matrices 
elas <- (A/L1) * S 
elasticity <- round(elas, 3)
