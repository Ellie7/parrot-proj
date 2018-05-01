####### YSA demography 'scraps'
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

