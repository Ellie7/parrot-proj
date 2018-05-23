# Yellow-shouldered amazon elasticity analyses 
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
source(file = "YSA life history data.R")

#inputting data
# mean hatchability and mean nestling survival with sample size and standard error 
total_summary

#manually creating data frame for now
#creating columns
stage <- c("1a", "1b", "1c", "2", "3")
# stage classes 
class <- c("egg", "nestling", "fledgling", "juvenile", "adult")
#stage durations 
di <- c((27/365), (59/365), (279/365), (24/12), 10) #27 days, 59 days, To age 12 months, Age 13-36 months, 
# Age 37 months+ (as 10 years -> needs updating and a reference) 
#pi (survival)
pi <- c((total_summary$mean_hatch[1]), (total_summary$mean_nestling_surv[1]), 0.71, 0.875, 0.875) 
# ^from life_table_data_master_csv, 0.71 from Salinas-Melgoza & Renton 2007, 0.875 from Rodriguez et al 2004 
# pi standard errors / SD 
piSD <- c((total_summary$se_hatch[1]), (total_summary$se_nestling_surv[1]), 0.2, 0.075, 0.075) # 
# ^using life_table_data_master_csv, 0.2 from Salinas-Melgoza & Renton 2007, 0.075s from Rodriguez et al 2004 
#reproductive output/fecundity 
f <- c(0, 0, 0, 0, 1.6) #half of 3.2 as sex ratio assumed 1:1 (Sams thesis)
#reproductive output/fecundity SEs
fSD <- c(0, 0, 0, 0, 0.1) # SE sams thesis (halved)
# creating dataframe by combining columns 
yellow_orig <- data_frame(stage, class, di, pi, piSD,  f, fSD)

#-----------------------------------------------------------------------------------------------------------------------------------

# FINAL MEAN DATA FRAME (yellow) for analysis 

#creating columns
stage <- c("1a", "1b", "1c", "2", "3")
# stage classes 
class <- c("egg", "nestling", "fledgling", "juvenile", "adult")
#stage durations 
di <- c((27/365), (59/365), (279/365), (24/12), 10) #27 days, 59 days, To age 12 months, Age 13-36 months, 
# Age 37 months+ (as 10 years -> needs updating and a reference) 
#pi (survival)
pi <- c((total_summary$mean_hatch[1]), (total_summary$mean_nestling_surv[1]), 0.71, 0.8515058, 0.8515058) 
# ^from life_table_data_master_csv, 0.71 from Salinas-Melgoza & Renton 2007, 0.8515s from Tamora's imputation 
# pi standard errors / SD 
piSD <- c((total_summary$se_hatch[1]), (total_summary$se_nestling_surv[1]), 0.2, 0.04722409, 0.04722409)
# ^ using life_table_data_master_csv, 0.2 from Salinas-Melgoza & Renton 2007, 0.1322s from Tamora's imputation  
#^ SD generated from the imputation (0.36363) was too high to generate beta values so you replaced it with a value slightly less 
# than the maximum allowed value (0.3555892).
#reproductive output/fecundity 
f <- c(0, 0, 0, 0, 1.6) #half of 3.2 as sex ratio assumed 1:1 (Sams thesis)
#reproductive output/fecundity SEs
fSD <- c(0, 0, 0, 0, 0.1) # SE sams thesis (halved)
# creating dataframe by combining columns 
yellow <- data_frame(stage, class, di, pi, piSD,  f, fSD)

#initial calculations using ysaFunc & ysameanFunc to create matrices 
A <- ysameanFunc(yellow) # for 'mean' matrix 
B <- ysaFunc(yellow) #for matrix drawn randomly from beta and lognormal distributed vital rates 

#eigen analysis 
eigs.A <- eigen(A)
eigs.A

#finding the first eigenvalue (finite rate of increase), lambda 
dom.pos <- which.max(eigs.A[["values"]])
L1mean <- Re(eigs.A[["values"]][dom.pos])
L1mean
lambda <- Re(eigs.A$values[1])
#=0.9329156

#finding r (log of lambda)
r <- log(L1mean)
r

#calculating the stable stage distribution 
w <- Re(eigs.A[["vectors"]][, dom.pos])
ssd <- w / sum(w)
stable <- ssd*100
stable <- round(stable, 2)
#0.207 0.670 0.115 0.007 0.000 0.000 0.002 
#calculating the reproductive value 
M <- eigen(t(A))
v <- Re(M$vectors[,which.max(Re(M$values))])
RV <- v / v[1]
RV 

#### or use  eigen.analysis from popbio 
eigen.analysis(A)


#------------------------------------------------------------------
#creating table of stable stage distribution & reproductive values 
stage_number <- c(1,2,3)
class_label <- c("Egg", "Juvenile", "Adult")
tab_5<- data.frame(stage_number, class_label, stable, RV)
colnames(tab_5) <- c("Stage number", "Stage Class", "Stable stage distribution (Dominant eigenvector)", "Reproductive values 
                     (left eigenvector)" )
tab_5
kable(tab_5, caption = "Table 1. Stable stage distribution (wJ) and reproductive values (v') for the yellow shouldered amazon 
      population matrix.")

#creating a theme for my graphs so that they are uniform in presentation 
mytheme <- theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + 
  theme(axis.title = element_text(size = 14))

#--------------------------------------------------------------------------------------------------------------------------------
# creating a figure which shows changes in rate of increase r resulting from simulated changes (10%) in fecundity and survival of 
#individual 
#life history stages  (similar to figure 1a & 1b in crouse 1987
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
yellow_pi3Inc <- mutate(yellow, pi = ifelse(stage == "3",  max(pi * 1.1, 0.99), pi *1)) #can't handle 1s again don't know why
#decreases 
AFd <- ysameanFunc(yellow_fecAdjust) 
eigs.AFd <- eigen(AFd)
A1ad <- ysameanFunc(yellow_pi1aAdj) 
eigs.A1ad <- eigen(A1ad)
A1bd <- ysameanFunc(yellow_pi1bAdj) 
eigs.A1bd <- eigen(A1bd)
A1cd <- ysameanFunc(yellow_pi1cAdj) 
eigs.A1cd <- eigen(A1cd)
A2d <- ysameanFunc(yellow_pi2Adj) 
eigs.A2d <- eigen(A2d)
A3d <- ysameanFunc(yellow_pi3Adj) 
eigs.A3d <- eigen(A3d)
#increases 
AFi <- ysameanFunc(yellow_fecIncrease) 
eigs.AFi <- eigen(AFi)
A1ai <- ysameanFunc(yellow_pi1aInc) 
eigs.A1ai <- eigen(A1ai)
A1bi <- ysameanFunc(yellow_pi1bInc) 
eigs.A1bi <- eigen(A1bi)
A1ci <- ysameanFunc(yellow_pi1cInc) 
eigs.A1ci <- eigen(A1ci)
A2i <- ysameanFunc(yellow_pi2Inc) 
eigs.A2i <- eigen(A2i)
A3i <- ysameanFunc(yellow_pi3Inc) 
eigs.A3i <- eigen(A3i)
#finding the first eigenvalue (finite rate of increase)
#decreases
dom.pos.Fd <- which.max(eigs.AFd[["values"]])
LFd <- Re(eigs.AFd[["values"]][dom.pos.Fd]) 
dom.pos.1ad <- which.max(eigs.A1ad[["values"]])
L1ad <- Re(eigs.A1ad[["values"]][dom.pos.1ad]) 
dom.pos.1bd <- which.max(eigs.A1bd[["values"]])
L1bd <- Re(eigs.A1bd[["values"]][dom.pos.1bd]) 
dom.pos.1cd <- which.max(eigs.A1cd[["values"]])
L1cd <- Re(eigs.A1cd[["values"]][dom.pos.1cd]) 
dom.pos.2d <- which.max(eigs.A2d[["values"]])
L2d <- Re(eigs.A2d[["values"]][dom.pos.2d]) 
dom.pos.3d <- which.max(eigs.A3d[["values"]])
L3d <- Re(eigs.A3d[["values"]][dom.pos.3d]) 
#increases 
dom.pos.Fi <- which.max(eigs.AFi[["values"]])
LFi <- Re(eigs.AFi[["values"]][dom.pos.Fi]) 
dom.pos.1ai <- which.max(eigs.A1ai[["values"]])
L1ai <- Re(eigs.A1ai[["values"]][dom.pos.1ai]) 
dom.pos.1bi <- which.max(eigs.A1bi[["values"]])
L1bi <- Re(eigs.A1bi[["values"]][dom.pos.1bi]) 
dom.pos.1ci <- which.max(eigs.A1ci[["values"]])
L1ci <- Re(eigs.A1ci[["values"]][dom.pos.1ci]) 
dom.pos.2i <- which.max(eigs.A2i[["values"]])
L2i <- Re(eigs.A2i[["values"]][dom.pos.2i]) 
dom.pos.3i <- which.max(eigs.A3i[["values"]])
L3i <- Re(eigs.A3i[["values"]][dom.pos.3i]) 

#plotting a figure -  Changes in rate of increase r resulting from simulated changes (10%) in fecundity and survival of individual 
#life history stages  (similar to figure 1a & 1b in crouse 1987)
lambdas<- c(LFd, L1ad, L1bd, L1cd, L2d, L3d)
rsa <- log(lambdas)
lambdas<- c(LFi, L1ai, L1bi, L1ci, L2i, L3i)
rsb <- log(lambdas)
#have a look
r
rsa
rsb
#data frame
r <- c(-(r-rsa[1]), -(r-rsa[2]), -(r-rsa[3]), -(r-rsa[4]), -(r-rsa[5]), -(r-rsa[6]), (rsb[1]-r), (rsb[2]-r),(rsb[3]-r), (rsb[4]-r), (rsb[5]-r), (rsb[6]-r))
change <- c("(a) Decrease 10%", "(a) Decrease 10%", "(a) Decrease 10%", "(a) Decrease 10%", "(a) Decrease 10%","(a) Decrease 10%",
            "(b) Increase 10%","(b) Increase 10%","(b) Increase 10%","(b) Increase 10%","(b) Increase 10%","(b) Increase 10%")
class <- c("Fecundity","Egg","Nestling","Fledgling","Juvenile","Adult","Fecundity","Egg","Nestling","Fledgling","Juvenile","Adult")
bar <- data_frame(class, change, r)

#ggplot 
figbar <- ggplot(bar, aes(x=class, y=r)) + geom_bar(stat = "identity") + facet_wrap(~change)
figbar <- figbar + mytheme
figbar <- figbar + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
figbar <- figbar + scale_x_discrete(limits=c("Fecundity","Egg","Nestling","Fledgling","Juvenile","Adult"))
figbar <- figbar + labs(x = "Stage class", y = "Change in r", size = 20)
figbar

#--------------------------------------------------------------------------------------------------------------------------------
# similar to crouse figure 3 
# creating a figure which shows the elasticity, or proportional sensitivity, of lambda to changes in fecundity F, survival 
# while remaining in the same stage P, and survival with growth, G. Because the elasticities of these matrix elements sum to 1, 
# they can be compared directly in terms of their contribution to the population growth rate 
#eigen analysis 
A <- ysameanFunc(yellow)
eigs.A <- eigen(A)
eigs.A
#finding the first eigenvalue (finite rate of increase)
dom.pos <- which.max(eigs.A[["values"]])
L1mean <- Re(eigs.A[["values"]][dom.pos])
L1mean

#sensitivity of projection matrices 
vw.s <- v %*% t (w) 
S <- vw.s/as.numeric(v %*% w)

#elasticity of projection matrices 
elasticity <- (A/L1mean) * S 


### plotting figure 3 - plot the proportional sensitivity to changes in F, P and G 
stage <- c("E", "J", "A")
#F <- c(elasticity[1, 1:3]) if including F1 & F2 on graph 
F <- c(NA, NA, elasticity[1,3])
#P <- c(elasticity[1,1], elasticity[2,2], elasticity[3,3]) if including P1 
P <- c(NA, elasticity[2,2], elasticity[3,3])
G <- c(elasticity[2,1], elasticity[3,2], NA)
sensitivities <- data.frame(stage, F, P, G)
sensitivity <- read.csv("cheat for now.csv") 
sens <- gather(sensitivities, vr, elasticity, F, P, G)
mutate(sens, "vr" = "Vital rate")
#change sensitivites data frame into the correct format 
fig <- ggplot(sens, aes(x = stage, y = elasticity, colour = vr, vr)) + geom_point(size = 6) + labs(x = "Stage", y = "Elasticity", size = 20)
fig <- fig + mytheme
fig <- fig + scale_x_discrete(limits=c("E","J","A"), labels=c("Egg", "Juvenile", "Adult"))
fig<- fig + ylim(0, max(elasticity))
fig <-fig + scale_colour_discrete(name="Vital rate",
                              breaks=c("F", "G", "P"),
                              labels=c("Reproductive Output", "Survival + growth", "Survival"))
fig <- fig + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())
fig



#--------------------------------------------------------------------------------------------------------------------------------
# similar to crouse figure 2 
# Looking at the effect of Age of First Reproduction on Intrinsic rate of Increase (r)
#-------------- Calculating changes in rate of increase r resulting from subtracting and adding 1 yr to the 
#calculations of P, and G, for each of the three immature stages.
yellow_di1 <- mutate(yellow, di = ifelse(stage == "2", di + 1, di * 1))
yellow_dii1 <- mutate(yellow_di1, di = ifelse(stage == "3", di - 1, di * 1)) # stage 2 as increasing age at which first breeding is = increasing number of years as juvenile
yellow_di2 <- mutate(yellow, di = ifelse(stage == "2", di + 2, di * 1)) # ^therefore juvenile = birds which have fledged the nest but have not began attempts at breeding
yellow_dii2 <- mutate(yellow_di2, di = ifelse(stage == "3", di - 2, di * 1)) 
yellow_di3 <- mutate(yellow, di = ifelse(stage == "2", di + 3, di * 1)) #dii1 = di is stage duration, i1 means increased by 1
yellow_dii3 <- mutate(yellow_di3, di = ifelse(stage == "3", di - 3, di * 1))


#applying my matrix function to the 3 new tables 
matinc1 <- ysameanFunc(yellow_dii1)
matinc2 <- ysameanFunc(yellow_dii2)
matinc3 <- ysameanFunc(yellow_dii3)

#eigenanalyses 
eigs.matinc1 <- eigen(matinc1)
eigs.matinc2 <- eigen(matinc2)
eigs.matinc3 <- eigen(matinc3)


#finding the first eigenvalue (finite rate of increase)
dom.pos.i1 <- which.max(eigs.matinc1[["values"]])
L1 <- Re(eigs.matinc1[["values"]][dom.pos.i1])
dom.pos.i2 <- which.max(eigs.matinc2[["values"]])
L2 <- Re(eigs.matinc2[["values"]][dom.pos.i2])
dom.pos.i3 <- which.max(eigs.matinc3[["values"]])
L3 <- Re(eigs.matinc3[["values"]][dom.pos.i3])

#finding r (log of lambda)
r1 <- log(L1)
r2 <- log(L2)
r3 <- log(L3)


#plotting the figure 
age <- c(3, 4, 5, 6)
lambdas<- c(L1mean, L1, L2, L3)
rs <- log(lambdas)
table_agerep <- data.frame(age, rs)
figure2 <- ggplot(table_agerep, aes(x = age, y = rs)) + geom_line(size=1) + geom_point(size=2)
figure2 <- figure2 + mytheme
figure2 <- figure2 + labs(x = "Age of First Reproduction (yr)", y = "Intrinsic rate of Increase (r)")  
figure2 <- figure2 + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  annotate("text", x=3.3, y=0.055, label = "base run") +
  theme(axis.line = element_line(colour = "black")) 
figure2

#-----------------------------------------------------------------------------------------------------------------------------------
# some popbio stuff 

image2(A) #mean matrix 

