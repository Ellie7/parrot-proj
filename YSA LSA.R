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

n0 <- c(1, 1, 1)
stoch.projection(mat1, n0, tmax = 50, nreps = 10, prob = NULL,
                 nmax = NULL, sumweight = rep(1, length(n0)), verbose=FALSE)

####################################### alternatively use map and eigen.analysis 
# load data (list of 10 matrices)
mat1 <- map(1:10, function(x) ysaFunc(yellow))
names(mat1) <- paste('M', 1:10, sep = '')
# use map and eigen.analysis
out <- map(mat1, eigen.analysis) 

# the eigen analysis outputs for each matrix
out

#----------------------------------------------------------------------------------------------------------------------------------
#sensitivities 
# need this to get the names onto the matrix at the end....
# only for old skool barplotting.
use.names.r<- rownames(out[[1]]$sensitivities)
use.names.c<- colnames(out[[1]]$sensitivities)

# create a data frame where columns are the matrix ID and rows are the elements
# of the sensitivity matrix

out2s <- map_df(out, function(x) {cbind(c(x$sensitivities))}) 

# use rowMeans to get the average sensitivity for each element
# and the se_fnc to get the se
# and re-create as a matrix

se_fnc <- function(x){sd(x)/sqrt(sum(!is.na(x)))}

# these are the mean and se matrices of sensitivity, in this case.
mean_sensitivity <- matrix(rowMeans(out2s),3,3,
                           dimnames = list(use.names.r,use.names.c))
se_sensitivity <- matrix(apply(out2s, 1, function(x) se_fnc(x)), 3,3)

# ggplot figure making
mat_element <- c(NA, "G1", NA, NA, "P2", "G2", "F3", NA, "P3")

g_plot_dfs <- data.frame(expand.grid(stage_A = use.names.r, stage_B = use.names.c),
                        mean_sens = c(mean_sensitivity), se_sens = c(se_sensitivity), mat_element)
glimpse(g_plot_dfs)

ggplot(g_plot_dfs, aes(x = stage_A, y = mean_sens, group = stage_B, fill = stage_B))+
  geom_col()+
  facet_wrap(~stage_B) + 
  geom_errorbar(aes(ymin = mean-se_sens, ymax = mean+se_sens)) 

#without NAs & not facet wrapped 
clean_s <- na.omit(g_plot_dfs) 
fig_s <- ggplot(clean, aes(x = mat_element, y = mean_sens, fill = stage_B)) + geom_col() + 
  scale_x_discrete(limits=c("G1", "P2", "G2", "P3", "F3")) +
  labs(x = "Matrix Element", y = "Mean Sensitivity") 
fig_s + scale_fill_discrete(name="Stage Class",
                        breaks=c("col1", "col2", "col3"),
                        labels=c("Egg", "Juvenile", "Adult")) + 
  geom_errorbar(aes(ymin = mean_sens-se_sens, ymax = mean_sens+se_sens))

#----------------------------------------------------------------------------------------------------------------------------------
# elasticities 
use.names.r<- rownames(out[[1]]$elasticities)
use.names.c<- colnames(out[[1]]$elasticities)

# create a data frame where columns are the matrix ID and rows are the elements
# of the sensitivity matrix

out2e <- map_df(out, function(x) {cbind(c(x$elasticities))}) 

# use rowMeans to get the average elasticity for each element
# and the se_fnc to get the se
# and re-create as a matrix

se_fnc <- function(x){sd(x)/sqrt(sum(!is.na(x)))}

# these are the mean and se matrices of elasticity
mean_elasticity <- matrix(rowMeans(out2e),3,3,
                 dimnames = list(use.names.r, use.names.c))
se_elasticity <- matrix(apply(out2e, 1, function(x) se_fnc(x)), 3,3)

# ggplot figure making
mat_element <- c(NA, "G1", NA, NA, "P2", "G2", "F3", NA, "P3")
g_plot_dfe <- data.frame(expand.grid(stage_A = use.names.r, stage_B = use.names.c),
                        mean_elas = c(mean_elasticity), se_elas = c(se_elasticity), mat_element)

ggplot(g_plot_df, aes(x = stage_A, y = mean_elas, group = stage_B, fill = stage_B))+
  geom_col()+
  facet_wrap(~stage_B)
+ 
  geom_errorbar(aes(ymin = mean_elas-se_elas, ymax = mean_elas+se_elas)) 

#without NAs
clean_e <- na.omit(g_plot_dfe) 

fig_e <- ggplot(clean_e, aes(x = mat_element, y = mean_elas, fill = stage_B)) + geom_col()+ 
  scale_x_discrete(limits=c("G1", "P2", "G2", "P3", "F3")) +
  labs(x = "Matrix Element", y = "Mean Elasticity") 

fig_e <- fig_e + scale_fill_discrete(name="Stage Class",
                          breaks=c("col1", "col2", "col3"),
                          labels=c("Egg", "Juvenile", "Adult")) + 
  geom_errorbar(aes(ymin = mean_elas-se_elas, ymax = mean_elas+se_elas)) 

#---------------------------------------------------------------------------------------------------------------------------------- 
# lambda
use.names <- rownames(out[[1]]$lambda1)

# create a data frame where columns are the matrix ID and rows are the elements
# of the sensitivity matrix

out2l <- map_df(out, function(x) {cbind(c(x$lambda1))}) 

# use rowMeans to get the average lambda
# and the se_fnc to get the se
# and re-create as a matrix

se_fnc <- function(x){sd(x)/sqrt(sum(!is.na(x)))}

# these are the mean and se matrices of lambda
mean_lambda <- c(matrix(rowMeans(out2e),dimnames = list(use.names,use.names)))
se_lambda <- matrix(apply(out2l, 1, function(x) se_fnc(x)))

# ggplot figure making
c(P1, )
g_plot_dfe <- data.frame(expand.grid(stage_A = use.names, stage_B = use.names),
                         mean_elas = c(mean_elasticity), se_elas = c(se_elasticity))

ggplot(g_plot_df, aes(x = stage_A, y = mean_elas, group = stage_B, fill = stage_B))+
  geom_col()+
  facet_wrap(~stage_B) + 
  geom_errorbar(aes(ymin = mean-se_elas, ymax = mean+se_elas)) 