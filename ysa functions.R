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

# Here is the data that will input into the function
yellow <- data_frame(stage, class, di, pi, piSD,  f, fSD)

# function 
ysaFunc <- function (dataSource) 
{ 
  #ps
  p1a<- betaval((dataSource$pi[1]), (dataSource$piSD[1]), fx=runif(1)) 
  p1b<- betaval((dataSource$pi[2]), (dataSource$piSD[2]), fx=runif(1)) 
  p1c<- betaval((dataSource$pi[3]), (dataSource$piSD[3]), fx=runif(1))
  p2 <- betaval((dataSource$pi[4]), (dataSource$piSD[4]), fx=runif(1))
  p3 <- betaval((dataSource$pi[5]), (dataSource$piSD[5]), fx=runif(1))
  
  #f
  f3 <- rnorm(1, mean = (dataSource$f[5]), sd = (dataSource$fSD[5])) #should 3.3 be divided by 2  
  
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
  matrix2[1,1] <- ((1 - (p1^(d1 - 0.99)))/(1 - (p1^d1)))*p1 #0.99 as doesn't like 1s 
  matrix2[2,2] <- ((1 - (p2^(d2 - 1)))/(1 - (p2^d2)))*p2
  matrix2[3,3] <- ((1 - (p3^(d3 - 1)))/(1 - (p3^d3)))*p3
  
  #add f
  matrix2[1,3] <- f3
  
  #add gs 
  matrix2[2,1] <- (p1^d1*(1 - p1))/(1 - p1^d1) 
  matrix2[3,2] <- (p2^d2*(1 - p2))/(1 - p2^d2) 
  
  return(matrix2)
} 

ysaFunc()
