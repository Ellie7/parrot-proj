# Yellow-shouldered amazon life stage simulation analysis 
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
source(file = "ysa functions.R")
source(file = "ysa demography.R") 

# making matrices using 'ysaFunc' function 
ysaFunc(yellow)
yellow
#make 10 matrices of random numbers.#for current code these aren't random and are drawn from ysa beta and normal distributed 
#vital rates
mat1 <- map(1:10, function(x) ysaFunc(yellow))

# use these matrices in mat 1 and get the eigen system for each....
mat2<-map(mat1, function(x) eigen(x)) 

#finding the first eigenvalue (finite rate of increase), lambda for each matrix 
mat3 <- map(mat2, function(x) {
  dom.pos <- which.max(x[["values"]])
  Re(x[["values"]][dom.pos])
})
#rs 
mat4 <- map(mat3, function(x) log(x)) 

##### alternatively use map and eigen.analysis 
analysis <- map(mat1, eigen.analysis) 

       
n0 <- c(1, 1, 1)
stoch.projection(mat1, n0, tmax = 50, nreps = 10, prob = NULL,
                        nmax = NULL, sumweight = rep(1, length(n0)), verbose=FALSE)
# load hudsonia data (list of 4 matrices)
data('hudsonia')
# use map and eigen.analysis
out<-map(hudsonia, eigen.analysis)

# the eigen analysis outputs for each matrix
out

# need this to get the names onto the matrix at the end....
# only for old skool barplotting.
# for hudsonia, its seed, seedlings, tiny etc....
use.names <- rownames(out[[1]]$sensitivities)

# create a data frame where columns are the matrix ID and rows are the elements
# of the sensitivity matrix

out2 <- map_df(out, function(x) {cbind(c(x$sensitivities))})

# use rowMeans to get the average sensitivity for each element
# and the se_fnc to get the se
# and re-create as a matrix

se_fnc <- function(x){sd(x)/sqrt(sum(!is.na(x)))}

# these are the mean and se matrices of sensitivity, in this case.
mean_sensitivity <- matrix(rowMeans(out2),6,6,
                           dimnames = list(use.names,use.names))
se_sensitivity <- matrix(apply(out2, 1, function(x) se_fnc(x)), 6,6)

# ggplot figure making
g_plot_df <- data.frame(expand.grid(stage_A = use.names, stage_B = use.names),
                        mean_sens = c(mean_sensitivity), se_sens = c(se_sensitivity))

ggplot(g_plot_df, aes(x = stage_A, y = mean_sens, group = stage_B, fill = stage_B))+
  geom_col()+
  facet_wrap(~stage_B)


    