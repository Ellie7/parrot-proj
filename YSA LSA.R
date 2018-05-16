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


#---------------------------------------------------------------------------------------------------------------------------------- 
# age at first reproduction plot with error bars 

#calculations of P, and G, for each of the three immature stages.
yellow_di1 <- mutate(yellow, di = ifelse(stage == "2", di + 1, di * 1))
yellow_dii1 <- mutate(yellow_di1, di = ifelse(stage == "3", di - 1, di * 1)) # stage 2 as increasing age at which first breeding is = increasing number of years as juvenile
yellow_di2 <- mutate(yellow, di = ifelse(stage == "2", di + 2, di * 1)) # ^therefore juvenile = birds which have fledged the nest but have not began attempts at breeding
yellow_dii2 <- mutate(yellow_di2, di = ifelse(stage == "3", di - 2, di * 1)) 
yellow_di3 <- mutate(yellow, di = ifelse(stage == "2", di + 3, di * 1)) #dii1 = di is stage duration, i1 means increased by 1
yellow_dii3 <- mutate(yellow_di3, di = ifelse(stage == "3", di - 3, di * 1))


#load data, creating groups of matrices with different ages at first breeding 
mat_inc1 <- map(1:10, function(x) ysaFunc(yellow_dii1))
mat_inc2 <- map(1:10, function(x) ysaFunc(yellow_dii2))
mat_inc3 <- map(1:10, function(x) ysaFunc(yellow_dii3))
names(mat_inc1) <- paste('M', 1:10, sep = '')
names(mat_inc2) <- paste('M', 1:10, sep = '')
names(mat_inc3) <- paste('M', 1:10, sep = '')

# use map and eigen.analysis
out_1 <- map(mat_inc1, eigen.analysis) 
out_2 <- map(mat_inc2, eigen.analysis)
out_3 <- map(mat_inc3, eigen.analysis)
# the eigen analysis outputs for each matrix
out_1
out_2
out_3 

out2_1 <- map_df(out_1, function(x) {cbind(c(x$lambda1))})
out2_2 <- map_df(out_2, function(x) {cbind(c(x$lambda1))})
out2_3 <- map_df(out_3, function(x) {cbind(c(x$lambda1))})


# these are the mean and se matrices of lambda
mean_lambda_1 <- c(matrix(rowMeans(out2_1)))
se_lambda_1 <- matrix(apply(out2_1, 1, function(x) se_fnc(x)))
mean_lambda_2 <- c(matrix(rowMeans(out2_2)))
se_lambda_2 <- matrix(apply(out2_2, 1, function(x) se_fnc(x)))
mean_lambda_3 <- c(matrix(rowMeans(out2_3)))
se_lambda_3 <- matrix(apply(out2_3, 1, function(x) se_fnc(x)))

# plotting graph
#making a dataframe to plot 

age <- c(3, 4, 5, 6)
lambdas<- c(mean_lambda, mean_lambda_1, mean_lambda_2, mean_lambda_3)
SE <- c(se_lambda, se_lambda_1, se_lambda_2, se_lambda_3)
rs <- log(lambdas)
rsSE <- log(SE)


# lambda vs repro graph 
table_agerep <- data.frame(age, lambdas, SE)
figure2 <- ggplot(table_agerep, aes(x = age, y = lambdas)) + geom_line(size=1) + geom_point(size=3)
figure2 <- figure2 + labs(x = "Age of First Reproduction (yr)", y = "Population Growth (lambda)")+ 
  geom_errorbar(aes(ymin = lambdas-SE, ymax = lambdas+SE), width = 0.1)  
figure2 <- figure2 + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  geom_vline(xintercept = 3, linetype = "dashed") +
  annotate("text", x=3.3, y=1.22, label = "base run") +
  theme(axis.line = element_line(colour = "black")) + theme(axis.title = element_text(size = 14)) 
figure2




