### figure 1 (a)
#changes in rate of increase r resulting from simulated changes in fecundity and survival of individual life history 
#caculating r determined in the baseline run of the matrix
table.3.50 <- read.csv("~/1 UNIVERSITY/Level 4/Project & Dissertation/Crouse 1987/table 3 surv fec altered .csv")
View(table.3.50)
#first recalculate eigenvalues for a 50% decrease in survivorship & 50% decrease in fecundity
A <- myFunc(table.3.50)
eigs.A <- eigen(A)
eigs.A 
#finding the first eigenvalue (finite rate of increase)
dom.pos <- which.max(eigs.A[["values"]])
L1N <- Re(eigs.A[["values"]][dom.pos]) #N for survivorship to distibuish from initial matrix 
L1N #=0.4046335  
lambda <- eigs.A[["values"]]
exp <- (Re(eigs.A[["values"]]))^2 #without doing ^2 for some reason get "Warning message: In log(lambdas) : NaNs produced"
rs <- log(sqrt(exp)) # sqrt to remove ^2 (see line above for explanation)
#------------------- Plotting Figure 1a
stage <- c("Eggs/Hatchlings", "Small Juveniles", "Large Juveniles", "Subadults", "Novice Breeders", "1st-yr Remigrants", "Mature Breeders")
changes <- data.frame(stage, rs)
ggplot(changes, aes(x = stage, y = rs)) + geom_bar(stat = "identity")+ labs(x = "Stage Class", y = "Intrinsic rate of Increase (r)") + geom_hline(aes(yintercept=-0.056399), linetype = "dashed", size = 0.75) + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ expand_limits(y=1)+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 
#^ if its not all on one line it errors its weird 
### Figure 1 (b)
table.3.increase <- read.csv("~/1 UNIVERSITY/Level 4/Project & Dissertation/Crouse 1987/table 3 surv fec increased .csv")
View(table.3.increase)
table.3.exp <- read.csv("~/1 UNIVERSITY/Level 4/Project & Dissertation/Crouse 1987/table 3 surv fec increased exp.csv")
View(table.3.exp)
#first recalculate eigenvalues for a 50% decrease in survivorship & 50% decrease in fecundity
B <- myFunc(table.3.exp)
eigs.B <- eigen(B)
eigs.B
#finding the first eigenvalue (finite rate of increase)
dom.pos <- which.max(eigs.B[["values"]])
L1N <- Re(eigs.B[["values"]][dom.pos]) #N for survivorship to distibuish from initial matrix 
L1N #=  
lambda <- eigs.B[["values"]]
exp <- (Re(eigs.B[["values"]]))^2 #without doing ^2 for some reason get "Warning message: In log(lambdas) : NaNs produced"
rs <- log(sqrt(exp)) # sqrt to remove ^2 (see line above for explanation)
#------------------- Plotting Figure 1 
stage <- c("Eggs/Hatchlings", "Small Juveniles", "Large Juveniles", "Subadults", "Novice Breeders", "1st-yr Remigrants", "Mature Breeders")
changes <- data.frame(stage, rs)
ggplot(changes, aes(x = stage, y = rs)) + geom_bar(stat = "identity")+ labs(x = "Stage Class", y = "Intrinsic rate of Increase (r)") + geom_hline(aes(yintercept=-0.056399), linetype = "dashed", size = 0.75) + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ expand_limits(y=1)+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 
#^ if its not all on one line it errors its weird