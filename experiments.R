#experiment 
sc1 <- rnorm(1, mean = (chicken$s[1]), sd = 0.0129)  
sc2 <- rnorm(1, mean = (chicken$s[5]), sd = (chicken$ssd[5])) 
sc3<- rnorm(1, mean = (chicken$s[6]), sd = (chicken$ssd[6]))
sc4 <- rnorm(1, mean = (chicken$s[7]), sd = (chicken$ssd[7]))
#m 
mc2 <- rnorm(1, mean = (chicken$m[5]), sd = (chicken$msd[5]))
mc3 <- rnorm(1, mean = (chicken$m[6]), sd = (chicken$msd[6]))
mc4 <- rnorm(1, mean = (chicken$m[7]), sd = (chicken$msd[7]))
mc5 <- rnorm(1, mean = (chicken$m[7]), sd = (chicken$msd[7]))

chickFunc <- function (chicken) 
{ sc1 <- rnorm(1, mean = (chicken$s[1]), sd = 0.0129)  
  sc2 <- rnorm(1, mean = (chicken$s[5]), sd = (chicken$ssd[5])) 
  sc3 <- rnorm(1, mean = (chicken$s[6]), sd = (chicken$ssd[6]))
  sc4 <- rnorm(1, mean = (chicken$s[7]), sd = (chicken$ssd[7]))
  #m 
  mc2 <- rnorm(1, mean = (chicken$m[5]), sd = (chicken$msd[5]))
  mc3 <- rnorm(1, mean = (chicken$m[6]), sd = (chicken$msd[6]))
  mc4 <- rnorm(1, mean = (chicken$m[7]), sd = (chicken$msd[7]))
  mc5 <- rnorm(1, mean = (chicken$m[7]), sd = (chicken$msd[7]))
  matrix2 <- matrix(0, nrow = 4, ncol = 4)
  #add sxmx
  matrix2[1,1] <- (sc1*mc2) 
  matrix2[1,2] <- (sc2*mc3)
  matrix2[1,3] <- (sc3*mc4)
  matrix2[1,4] <- (sc4*mc5)
  #add s(mxc$m[7]))
  matrix2[2,1] <- sc1
  matrix2[3,2] <- sc2
  matrix2[4,3] <- sc3
  matrix2[4,4] <- sc4
  return(matrix2)
}

chickFunc(chicken) #chickFunc returns a different matrix each time 
