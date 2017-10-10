#Crouse et al. 1987 
#A STAGE-BASED POPULATION MODEL FOR LOGGERHEAD SEA TURTLES AND IMPLICATIONS FOR CONSERVATION
#recreating table 4
A <- matrix(c(0, 0, 0, 0, 127, 4, 80, 0.6747, 0.7370, 0, 0,0, 0, 0, 0, 0.0486, 0.6610, 0, 0, 
              0, 0, 0, 0, 0.0147, 0.6907, 0, 0, 0, 0, 0, 0, 0.0518, 0, 0, 0, 0, 0, 0, 0, 0.8091, 
              0, 0, 0, 0, 0, 0, 0, 0.8091, 0.8089), nr=7, byrow = TRUE) 
#--------------- population projections 
#stage structure growth (multiple steps)
N0 <- 
years <- 
N.projections <- matrix(0, nrow = nrow(A), ncol = years + 1)
N.projections[,1] <- N0 
for (i in 1:years) N.projections[, i + 1] <- A %*% N.projections[, i]

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
for (i in 1:t) R.t[i] <- {
  Nt1 <- A %*% Nt
  R <- sum(Nt1)/sum(Nt)
  R
}
#compare result to the point estimate (L1)

#calculating the stable stage distribution 
w <- Re(eigs.A[["vectors"]][, dom.pos])
ssd <- w / sum(w)
round(ssd,3)
#0.207 0.670 0.115 0.007 0.000 0.000 0.002 
#calculating the reproductive value 
M <- eigen(t(A))
v <- Re(M$vectors[,which.max(Re(M$values))])
RV <- v / v[1]
RV 
#1.000000   1.400668   5.995523 115.844511 568.780852 507.373040 587.669314
#---------------- sensitivity analyses 