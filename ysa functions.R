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

ysaFunc <- function (dataSource) 
{ 
  #ps
  p1a<- betaval((pi[1]), (piSD[1]), fx=runif(1)) 
  p1b<- betaval((pi[2]), (piSD[2]), fx=runif(1)) 
  p1c<- betaval((pi[3]), (piSD[3]), fx=runif(1))
  p2 <- betaval((pi[4]), (piSD[4]), fx=runif(1))
  p3 <- betaval((pi[5]), (piSD[5]), fx=runif(1))
  #f
  f3 <- rnorm(1, mean = (f[5]), sd = (fSD[5])) #should 3.3 be divided by 2  
  # Pi <- ((1 - (pi^(di - 1)))/(1 - (pi^di)))*pi ------- equation for Pi's
  # Gi <- (pi^di*(1 - pi))/(1 - pi^di)           ------- equation for Gi's
  matrix2 <- matrix(0, nrow = 3, ncol = 3)
  d1 <- di[1] + di[2] + di[3]
  d2 <- di[4]
  d3 <- di[5]
  p1 <- (p1a*p1b*p1c) # this stage as the survival is from the multiplication of  p1a, p1b and p1c
  #add ps 
  matrix2[1,1] <- ((1 - (p1^(d1 - 0.99)))/(1 - (p1^d1)))*p1 #0.99 as doesn't lie 1s 
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