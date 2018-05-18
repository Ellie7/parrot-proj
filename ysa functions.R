#Yellow-shouldered amazon functions 
library(popbio)
library(tidyverse) 

# Here is the function
# pi, piSD etc are inside the data frame.  This needs to be specified this
# for example: dataSource$pi[1] is where pi[1] is located


#--------------------------------------------------------------------------------------------------------------------------------
# ysaFunc creates a matrix from randomly drawn data, with p values taken from a beta distribution and f values
# drawn from a lognormal distribution 

ysaFunc <- function (dataSource) 
{ 
  #ps
  p1a<- betaval((dataSource$pi[1]), (dataSource$piSD[1]), fx=runif(1)) 
  p1b<- betaval((dataSource$pi[2]), (dataSource$piSD[2]), fx=runif(1)) 
  p1c<- betaval((dataSource$pi[3]), (dataSource$piSD[3]), fx=runif(1))
  p2 <- betaval((dataSource$pi[4]), (dataSource$piSD[4]), fx=runif(1))
  p3 <- betaval((dataSource$pi[5]), (dataSource$piSD[5]), fx=runif(1))
  
  #f
  f3 <- rnorm(1, mean = (dataSource$f[5]), sd = (dataSource$fSD[5])) 
  
  # Pi <- ((1 - (pi^(di - 1)))/(1 - (pi^di)))*pi ------- equation for Pi's
  # Gi <- (pi^di*(1 - pi))/(1 - pi^di)           ------- equation for Gi's
  
  #d
  d1 <- dataSource$di[1] + dataSource$di[2] + dataSource$di[3]
  d2 <- dataSource$di[4]
  d3 <- dataSource$di[5]
  
  # this uses p1's defined above
  p1 <- (p1a*p1b*p1c) # this stage as the survival is from the multiplication of  p1a, p1b and p1c
  #add ps 
  
  # construct the matrix using defined parameters above
  matrix2 <- matrix(0, nrow = 3, ncol = 3)
  dimnames(matrix2) <- list(rownames(matrix2, do.NULL = FALSE, prefix = "row"),
                            colnames(matrix2, do.NULL = FALSE, prefix = "col"))
  matrix2[1,1] <- ((1 - (p1^(d1 - 1)))/(1 - (p1^d1)))*p1 
  matrix2[2,2] <- ((1 - (p2^(d2 - 1)))/(1 - (p2^d2)))*p2
  matrix2[3,3] <- ((1 - (p3^(d3 - 1)))/(1 - (p3^d3)))*p3
  
  #add f
  matrix2[1,3] <- f3
  
  #add gs 
  matrix2[2,1] <- (p1^d1*(1 - p1))/(1 - p1^d1) 
  matrix2[3,2] <- (p2^d2*(1 - p2))/(1 - p2^d2) 
  
  return(matrix2)
} 

#--------------------------------------------------------------------------------------------------------------------------------
# ysameanFunc creates a matrix based on means, ie the raw values from the data source rather than drawing from distributions

ysameanFunc <- function (dataSource) 
{ 
  #ps
  p1a<- dataSource$pi[1]
  p1b<- dataSource$pi[2]
  p1c<- dataSource$pi[3]
  p2 <- dataSource$pi[4]
  p3 <- dataSource$pi[5]
  
  #f
  f3 <- dataSource$f[5] #should 3.3 be divided by 2  
  
  # Pi <- ((1 - (pi^(di - 1)))/(1 - (pi^di)))*pi ------- equation for Pi's
  # Gi <- (pi^di*(1 - pi))/(1 - pi^di)           ------- equation for Gi's
  
  #d
  d1 <- dataSource$di[1] + dataSource$di[2] + dataSource$di[3]
  d2 <- dataSource$di[4]
  d3 <- dataSource$di[5]
  
  # this uses p1's defined above
  p1 <- (p1a*p1b*p1c) # this stage as the survival is from the multiplication of  p1a, p1b and p1c
  #add ps 
  
  # construct the matrix using defined parameters above
  matrix2 <- matrix(0, nrow = 3, ncol = 3)
  matrix2[1,1] <- ((1 - (p1^(d1 - 1)))/(1 - (p1^d1)))*p1  
  matrix2[2,2] <- ((1 - (p2^(d2 - 1)))/(1 - (p2^d2)))*p2
  matrix2[3,3] <- ((1 - (p3^(d3 - 1)))/(1 - (p3^d3)))*p3
  
  #add f
  matrix2[1,3] <- f3
  
  #add gs 
  matrix2[2,1] <- (p1^d1*(1 - p1))/(1 - p1^d1) 
  matrix2[3,2] <- (p2^d2*(1 - p2))/(1 - p2^d2) 
  
  return(matrix2)
} 

#--------------------------------------------------------------------------------------------------------------------------------
# ysaFuncDD creates a function to calculate a matrix model with density-dependent fecundity - need to pass in the current 
# population vector n and threshold for density-dependent effects. Also specify whether to make it stochastic or just use mean 
# values.
ysaFuncDD <- function (dataSource, n, threshold, stochastic = FALSE)
{
  if (stochastic) {
    #ps
    p1a<- betaval((dataSource$pi[1]), (dataSource$piSD[1]), fx=runif(1)) 
    p1b<- betaval((dataSource$pi[2]), (dataSource$piSD[2]), fx=runif(1)) 
    p1c<- betaval((dataSource$pi[3]), (dataSource$piSD[3]), fx=runif(1))
    p2 <- betaval((dataSource$pi[4]), (dataSource$piSD[4]), fx=runif(1))
    p3 <- betaval((dataSource$pi[5]), (dataSource$piSD[5]), fx=runif(1))
    
    # F
    # N > M  -> (M*F)/N
    # N <= M ->  F
    f3 <- rnorm(1, mean = (dataSource$f[5]), sd = (dataSource$fSD[5])) 
  } else {
    #ps
    p1a<- dataSource$pi[1]
    p1b<- dataSource$pi[2]
    p1c<- dataSource$pi[3]
    p2 <- dataSource$pi[4]
    p3 <- dataSource$pi[5]
    
    #f
    f3 <- dataSource$f[5]
  }
  
  # Pi <- ((1 - (pi^(di - 1)))/(1 - (pi^di)))*pi ------- equation for Pi's
  # Gi <- (pi^di*(1 - pi))/(1 - pi^di)           ------- equation for Gi's
  
  #d
  d1 <- dataSource$di[1] + dataSource$di[2] + dataSource$di[3]
  d2 <- dataSource$di[4]
  d3 <- dataSource$di[5]
  
  # this uses p1's defined above
  p1 <- (p1a*p1b*p1c) # this stage as the survival is from the multiplication of  p1a, p1b and p1c
  #add ps 
  
  # construct the matrix using defined parameters above
  matrix2 <- matrix(0, nrow = 3, ncol = 3)
  matrix2[1,1] <- ((1 - (p1^(d1 - 1)))/(1 - (p1^d1)))*p1  
  matrix2[2,2] <- ((1 - (p2^(d2 - 1)))/(1 - (p2^d2)))*p2
  matrix2[3,3] <- ((1 - (p3^(d3 - 1)))/(1 - (p3^d3)))*p3
  
  #add f - including density dependence based on number of breeders
  if (n[3] > threshold) {
    matrix2[1,3] <- f3*threshold/n[3]
  } else {
    matrix2[1,3] <- f3
  }
  
  #add gs 
  matrix2[2,1] <- (p1^d1*(1 - p1))/(1 - p1^d1) 
  matrix2[3,2] <- (p2^d2*(1 - p2))/(1 - p2^d2) 
  
  return(matrix2)
}
