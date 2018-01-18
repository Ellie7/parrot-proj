#code cut from main crouse because think its duplicated
#plotting figure 1a
stage_class <- c("Fecundity", "Eggs/Hatchlings", "Small Juveniles", "Large Juveniles", "Subadults", "Novice Breeders", "1st-yr Remigrants", "Mature Breeders")
lambdas<- c(LFd, L1d, L2d, L3d, L4d, L5d, L6d, L7d)
rs <- log(lambdas)
table_decrease <- data.frame(stage_class, rs)
graph <- ggplot(table_decrease, aes(x = stage_class, y = rs)) + geom_bar(stat = "identity")
graph1 <- graph + labs(x = "Stage Class", y = "Intrinsic rate of Increase (r)") + geom_hline(aes(yintercept=-0.056399), linetype = "dashed", size = 0.75)  
graph2 <- graph1 + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ expand_limits(y = -0.3, x = 0.4)
figure1a <- graph2 + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + scale_x_discrete(limits = c("Fecundity", "Eggs/Hatchlings", "Small Juveniles", "Large Juveniles", "Subadults", "Novice Breeders", "1st-yr Remigrants", "Mature Breeders"))
figure1a 
stage_class <- c("Fecundity", "Eggs/Hatchlings", "Small Juveniles", "Large Juveniles", "Subadults", "Novice Breeders", "1st-yr Remigrants", "Mature Breeders")
lambdas<- c(LFd, L1d, L2d, L3d, L4d, L5d, L6d, L7d)
rs <- log(lambdas)
table_decrease <- data.frame(stage_class, rs)
graph <- ggplot(table_decrease, aes(x = stage_class, y = rs)) + geom_bar(stat = "identity")
graph1 <- graph + labs(x = "Stage Class", y = "Intrinsic rate of Increase (r)") + geom_hline(aes(yintercept=-0.056399), linetype = "dashed", size = 0.75)  
graph2 <- graph1 + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ expand_limits(y = -0.3, x = 0.4)
figure1a <- graph2 + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + scale_x_discrete(limits = c("Fecundity", "Eggs/Hatchlings", "Small Juveniles", "Large Juveniles", "Subadults", "Novice Breeders", "1st-yr Remigrants", "Mature Breeders"))
figure1a 