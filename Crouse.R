#Crouse et al. 1987 
#A STAGE-BASED POPULATION MODEL FOR LOGGERHEAD SEA TURTLES AND IMPLICATIONS FOR CONSERVATION
install.packages("dplyr")
install.packages("ggplot2")
library(dplyr)
library(ggplot2)
##making the matrix: 
mat1 <- matrix(0, nrow = 7, ncol = 7, byrow = T)
#fecundity 
fecs <- select(table.3, fecundity)
mat1[1, 5] <- 127
mat1[1, 6] <- 4
mat1[1, 7] <- 80
#Gs
pi <- select(table.3, annual_survivorship)
di <-select(table.3, stage_duration)
Gi <- (pi^di*(1-pi))/(1-(pi^di)) 
round(Gi, 4)
#Ps
Pi <- (1-(pi^di))/(1-(pi^di))*pi
mat1[2:7, 1:6] <-Pi
#recreating table 4
A <- matrix(c(0, 0, 0, 0, 127, 4, 80, 0.6747, 0.7370, 0, 0,0, 0, 0, 0, 0.0486, 0.6610, 0, 0, 
              0, 0, 0, 0, 0.0147, 0.6907, 0, 0, 0, 0, 0, 0, 0.0518, 0, 0, 0, 0, 0, 0, 0, 0.8091, 
              0, 0, 0, 0, 0, 0, 0, 0.8091, 0.8089), nr=7, byrow = TRUE) 
#--------------- population projections 
#stage structure growth (multiple steps)
N0 <- 
years <- 10
N.projections <- matrix(0, nrow = nrow(A), ncol = years + 1) #^not happy yet
#annual growth rate
N.totals <- apply(N.projections, 2, sum)
Rs <- N.totals[-1]/N.totals[-(years + 1)]
plot(0:(years - 1), Rs, type = "b", xlab = "Year", ylab = "R") #^not happy yet
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
for (i in 1:t) R.t[i] <- {
  Nt1 <- A %*% Nt
  R <- sum(Nt1)/sum(Nt)
  R
} #^not happy until Nt is assigned 
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
#---------------------------------- sensitivity analyses 
#sensitivity of projection matrices 
vw.s <- v %*% t (w) 
(S <- vw.s/as.numeric(v %*% w))
#elasticity of projection matrices 
elas <- (A/L1) * S 
elasticity <- round(elas, 3)
