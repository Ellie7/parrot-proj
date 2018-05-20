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
fig_s <- ggplot(clean_s, aes(x = mat_element, y = mean_sens, fill = stage_B)) + geom_col() + 
  scale_x_discrete(limits=c("G1", "P2", "G2", "P3", "F3")) +
  labs(x = "Matrix Element", y = "Mean Sensitivity") 
fig_s <- fig_s + scale_fill_discrete(name="Stage Class",
                        breaks=c("col1", "col2", "col3"),
                        labels=c("Egg", "Juvenile", "Adult")) + 
  geom_errorbar(aes(ymin = mean_sens-se_sens, ymax = mean_sens+se_sens), width = 0.3)
fig_s + mytheme

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
  geom_errorbar(aes(ymin = mean_elas-se_elas, ymax = mean_elas+se_elas), width = 0.4) 
fig_e <- fig_e + mytheme
fig_e

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
mean_lambda <- c(matrix(rowMeans(out2l),dimnames = list(use.names,use.names)))

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
mat_inc1 <- map(1:100, function(x) ysaFunc(yellow_dii1))
mat_inc2 <- map(1:100, function(x) ysaFunc(yellow_dii2))
mat_inc3 <- map(1:100, function(x) ysaFunc(yellow_dii3))
names(mat_inc1) <- paste('M', 1:100, sep = '')
names(mat_inc2) <- paste('M', 1:100, sep = '')
names(mat_inc3) <- paste('M', 1:100, sep = '')

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
figure2 <- figure2 + mytheme
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme(axis.line = element_line(colour = "black")) + theme(axis.title = element_text(size = 14)) 
figure_2 <- figure2 + annotate("text", x=3.3, y=1.055, label = "base run") 
figure2

#---------------------------------------------------------------------------------------------------------------------------------- 
# sensitivity to pi i.e. the underlying vital rates

#decresing fecundity and survival by 10%
yellow_fecAdjust<-mutate(yellow, f = 0.9*f)
yellow_pi1aAdj <- mutate(yellow, pi = ifelse(stage == "1a", pi * 0.9, pi * 1)) 
yellow_pi1bAdj <- mutate(yellow, pi = ifelse(stage == "1b", pi * 0.9, pi * 1))
yellow_pi1cAdj <- mutate(yellow, pi = ifelse(stage == "1c", pi * 0.9, pi * 1))
yellow_pi2Adj <- mutate(yellow, pi = ifelse(stage == "2", pi * 0.9, pi * 1))
yellow_pi3Adj <- mutate(yellow, pi = ifelse(stage == "3", pi * 0.9, pi * 1))
#increase in fecundity or an increase in survivorship by 10%
yellow_fecIncrease<-mutate(yellow, f = 1.1*f)
yellow_pi1aInc <- mutate(yellow, pi = ifelse(stage == "1a", pi * 1.1, pi *1))
yellow_pi1bInc <- mutate(yellow, pi = ifelse(stage == "1b", pi * 1.1, pi *1))
yellow_pi1cInc <- mutate(yellow, pi = ifelse(stage == "1c", pi * 1.1, pi *1))
yellow_pi2Inc <- mutate(yellow, pi = ifelse(stage == "2", max(pi * 1.1, 0.99), pi *1)) #can't handle 1s again don't know why
yellow_pi3Inc <- mutate(yellow, pi = ifelse(stage == "3",  max(pi * 1.1, 0.99), pi *1))


#load data, creating groups of matrices with different ages at first breeding 
mat_decfds <- map(1:100, function(x) ysaFunc(yellow_fecAdjust))
mat_dec1ads <- map(1:100, function(x) ysaFunc(yellow_pi1aAdj))
mat_dec1bds <- map(1:100, function(x) ysaFunc(yellow_pi1bAdj))
mat_dec1cds <- map(1:100, function(x) ysaFunc(yellow_pi1cAdj))
mat_dec2ds <- map(1:100, function(x) ysaFunc(yellow_pi2Adj))
mat_dec3ds <- map(1:100, function(x) ysaFunc(yellow_pi3Adj))

mat_incfds <- map(1:100, function(x) ysaFunc(yellow_fecIncrease))
mat_inc1ads <- map(1:100, function(x) ysaFunc(yellow_pi1aInc))
mat_inc1bds <- map(1:100, function(x) ysaFunc(yellow_pi1bInc))
mat_inc1cds <- map(1:100, function(x) ysaFunc(yellow_pi1cInc))
mat_inc2ds <- map(1:100, function(x) ysaFunc(yellow_pi2Inc))
mat_inc3ds <- map(1:100, function(x) ysaFunc(yellow_pi3Inc))


names(mat_decfds) <- paste('M', 1:100, sep = '')
names(mat_dec1ads) <- paste('M', 1:100, sep = '')
names(mat_dec1bds) <- paste('M', 1:100, sep = '')
names(mat_dec1cds) <- paste('M', 1:100, sep = '')
names(mat_dec2ds) <- paste('M', 1:100, sep = '')
names(mat_dec3ds) <- paste('M', 1:100, sep = '')
names(mat_incfds) <- paste('M', 1:100, sep = '')
names(mat_inc1ads) <- paste('M', 1:100, sep = '')
names(mat_inc1bds) <- paste('M', 1:100, sep = '')
names(mat_inc1cds) <- paste('M', 1:100, sep = '')
names(mat_inc2ds) <- paste('M', 1:100, sep = '')
names(mat_inc3ds) <- paste('M', 1:100, sep = '')

# use map and eigen.analysis
out_df <- map(mat_decfds, eigen.analysis)    
out_d1a <- map(mat_dec1ads, eigen.analysis)
out_d1b <- map(mat_dec1bds, eigen.analysis)
out_d1c <-  map(mat_dec1cds, eigen.analysis)
out_d2 <-  map(mat_dec2ds, eigen.analysis)
out_d3 <-  map(mat_dec3ds, eigen.analysis)
out_if <- map(mat_incfds, eigen.analysis)    
out_i1a <- map(mat_inc1ads, eigen.analysis)
out_i1b <- map(mat_inc1bds, eigen.analysis)
out_i1c <-  map(mat_inc1cds, eigen.analysis)
out_i2 <-  map(mat_inc2ds, eigen.analysis)
out_i3 <-  map(mat_inc3ds, eigen.analysis)

# extracting senstivities from each matrix 
#decreases 
use.names.rdf<- rownames(out_df[[1]]$lambda1)
use.names.cdf<- colnames(out_df[[1]]$lambda1)
use.names.rd1a<- rownames(out_d1a[[1]]$lambda1)
use.names.cd1a<- colnames(out_d1a[[1]]$lambda1)
use.names.rd1b<- rownames(out_d1b[[1]]$lambda1)
use.names.cd1b<- colnames(out_d1b[[1]]$lambda1)
use.names.rd1c<- rownames(out_d1c[[1]]$lambda1)
use.names.cd1c<- colnames(out_d1c[[1]]$lambda1)
use.names.rd2<- rownames(out_d2[[1]]$lambda1)
use.names.cd2<- colnames(out_d2[[1]]$lambda1)
use.names.rd3<- rownames(out_d3[[1]]$lambda1)
use.names.cd3<- colnames(out_d3[[1]]$lambda1)
#increases
use.names.rif<- rownames(out_if[[1]]$lambda1)
use.names.cif<- colnames(out_if[[1]]$lambda1)
use.names.ri1a<- rownames(out_i1a[[1]]$lambda1)
use.names.ci1a<- colnames(out_i1a[[1]]$lambda1)
use.names.ri1b<- rownames(out_i1b[[1]]$lambda1)
use.names.ci1b<- colnames(out_i1b[[1]]$lambda1)
use.names.ri1c<- rownames(out_i1c[[1]]$lambda1)
use.names.ci1c<- colnames(out_i1c[[1]]$lambda1)
use.names.ri2<- rownames(out_i2[[1]]$lambda1)
use.names.ci2<- colnames(out_i2[[1]]$lambda1)
use.names.ri3<- rownames(out_i3[[1]]$lambda1)
use.names.ci3<- colnames(out_i3[[1]]$lambda1)
# create a data frame where columns are the matrix ID and rows are the elements
# of the sensitivity matrix

out2_df <- map_df(out_df, function(x) {cbind(c(x$lambda1))})    
out2_d1a <- map_df(out_d1a, function(x) {cbind(c(x$lambda1))}) 
out2_d1b <- map_df(out_d1b, function(x) {cbind(c(x$lambda1))})
out2_d1c <- map_df(out_d1c, function(x) {cbind(c(x$lambda1))})
out2_d2 <- map_df(out_d2, function(x) {cbind(c(x$lambda1))})
out2_d3 <- map_df(out_d3, function(x) {cbind(c(x$lambda1))})
out2_if <- map_df(out_if, function(x) {cbind(c(x$lambda1))})
out2_i1a <- map_df(out_i1a, function(x) {cbind(c(x$lambda1))})
out2_i1b <- map_df(out_i1b, function(x) {cbind(c(x$lambda1))})
out2_i1c <- map_df(out_i1c, function(x) {cbind(c(x$lambda1))})
out2_i2 <-  map_df(out_i2, function(x) {cbind(c(x$lambda1))})
out2_i3 <-  map_df(out_i3, function(x) {cbind(c(x$lambda1))})

# use rowMeans to get the average sensitivity for each element
# and the se_fnc to get the se
# and re-create as a matrix

se_fnc <- function(x){sd(x)/sqrt(sum(!is.na(x)))}

# these are the mean and se matrices of sensitivity, in this case.
#decreases 
df_mean_lambda <- c(matrix(rowMeans(out2_df)))
df_se_lambda <- matrix(apply(out2_df, 1, function(x) se_fnc(x)))
d1a_mean_lambda <- c(matrix(rowMeans(out2_d1a)))
d1a_se_lambda <- matrix(apply(out2_d1a, 1, function(x) se_fnc(x)))
d1b_mean_lambda <- c(matrix(rowMeans(out2_d1b)))
d1b_se_lambda <- matrix(apply(out2_d1b, 1, function(x) se_fnc(x)))
d1c_mean_lambda <- c(matrix(rowMeans(out2_d1c)))
d1c_se_lambda <- matrix(apply(out2_d1c, 1, function(x) se_fnc(x)))
d2_mean_lambda <- c(matrix(rowMeans(out2_d2)))
d2_se_lambda <- matrix(apply(out2_d2, 1, function(x) se_fnc(x)))
d3_mean_lambda <- c(matrix(rowMeans(out2_d3)))
d3_se_lambda <- matrix(apply(out2_d3, 1, function(x) se_fnc(x)))
#increases
if_mean_lambda <- c(matrix(rowMeans(out2_if)))
if_se_lambda <- matrix(apply(out2_if, 1, function(x) se_fnc(x)))
i1a_mean_lambda <- c(matrix(rowMeans(out2_i1a)))
i1a_se_lambda <- matrix(apply(out2_i1a, 1, function(x) se_fnc(x)))
i1b_mean_lambda <- c(matrix(rowMeans(out2_i1b)))
i1b_se_lambda <- matrix(apply(out2_i1b, 1, function(x) se_fnc(x)))
i1c_mean_lambda <- c(matrix(rowMeans(out2_i1c)))
i1c_se_lambda <- matrix(apply(out2_i1c, 1, function(x) se_fnc(x)))
i2_mean_lambda <- c(matrix(rowMeans(out2_i2)))
i2_se_lambda <- matrix(apply(out2_i2, 1, function(x) se_fnc(x)))
i3_mean_lambda <- c(matrix(rowMeans(out2_i3)))
i3_se_lambda <- matrix(apply(out2_i3, 1, function(x) se_fnc(x)))

# ggplot figure making
change_means<- c(-(mean_lambda-df_mean_lambda), -(mean_lambda-d1a_mean_lambda), -(mean_lambda-d1b_mean_lambda),  
                 -(mean_lambda-d1c_mean_lambda), -(mean_lambda-d2_mean_lambda), -(mean_lambda-d3_mean_lambda),
                 (if_mean_lambda - mean_lambda), (i1a_mean_lambda- mean_lambda), (i1b_mean_lambda - mean_lambda),  
                 (i1c_mean_lambda - mean_lambda), (i2_mean_lambda - mean_lambda), (i3_mean_lambda - mean_lambda))



change_se <- c(df_se_lambda, d1a_se_lambda, d1b_se_lambda,  d1c_se_lambda, d2_se_lambda, d3_se_lambda, 
               if_se_lambda, i1a_se_lambda, i1b_se_lambda,  i1c_se_lambda, i2_se_lambda, i3_se_lambda)

change <- c("decrease","decrease","decrease","decrease","decrease","decrease",
            "increase", "increase","increase","increase","increase","increase")
stage <- c("Fecundity", "Egg", "Nestling", "Fledgling", "Juvenile", "Adult",
           "Fecundity", "Egg", "Nestling", "Fledgling", "Juvenile", "Adult")

g_plot_under <- data_frame(stage,change, change_means, change_se)
lline <- (1-mean_lambda)
#ggplot 
figbar_LSA <- ggplot(g_plot_under, aes(x=stage, y=change_means)) + geom_bar(stat = "identity") + facet_wrap(~change)
figbar_LSA <- figbar_LSA + mytheme
figbar_LSA<- figbar_LSA + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
figbar_LSA <- figbar_LSA + scale_x_discrete(limits=c("Fecundity","Egg","Nestling","Fledgling","Juvenile","Adult"))
figbar_LSA <- figbar_LSA + labs(x = "Stage class", y = "Change in lambda", size = 20)
figbar_LSA <- figbar_LSA + geom_errorbar(aes(ymin = change_means-change_se, ymax = change_means+change_se), width = 0.2)  
figbar_LSA + geom_hline(yintercept = lline, linetype = "dashed") + annotate("text", x=3.3, y=-0.058, label = "lambda = 1") 
