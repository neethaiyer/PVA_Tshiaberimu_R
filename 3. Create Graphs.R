## Set the working directory:
workingDir <- "/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/"
setwd(workingDir)
workingDir_Results <- "/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_extn_results/"
##workingDir <- "~/Documents/git repositories/PVA_Tshiaberimu_R/"

## Let's pick our colors for MTN and WLG and include a transparency factor; note that these are only used for the bargraph below:
coral <- rgb(250, 128, 114, alpha=150, maxColorValue = 255) ## for WLG
azure4 <- rgb(131, 139, 139, alpha=150, maxColorValue = 255) ## for MTN
colfun <- colorRampPalette(c("azure4", "coral"))
colfun(4)

## A. Color for MTN with 3% growth rate
colMTN3 <- "#838B8B"
## B. Color for MTN with 2% growth rate
colMTN2 <- "#AC8777"
## C. Color for MTN with 1% growth rate
colMTN1 <- "#D58363"
## D. Color for WLG 
colWLG <- "#FF7F50"

##plot(rep(1,4), col=c(colMTN3,colMTN2,colMTN1,colWLG), pch=19 ,cex=3)

## Read all csv files:
dat <- read.csv(paste0(workingDir, "Gorilla_LifeTables.csv")) ## life tables
## probability of extinctions based on LM projections
wlg_lm <- read.csv(paste0(workingDir_Results,"extn_lm_WLG.csv"))
mtn_3per_lm <- read.csv(paste0(workingDir_Results,"extn_lm_MTN_3%.csv"))
mtn_2per_lm <- read.csv(paste0(workingDir_Results,"extn_lm_MTN_2%.csv"))
mtn_1per_lm <- read.csv(paste0(workingDir_Results,"extn_lm_MTN_1%.csv"))
wlg_ibm <- read.csv(paste0(workingDir_Results,"extn_ibm_WLG.csv"))
mtn_3per_ibm <- read.csv(paste0(workingDir_Results,"extn_ibm_MTN_3%.csv"))
mtn_2per_ibm <- read.csv(paste0(workingDir_Results,"extn_ibm_MTN_2%.csv"))
mtn_1per_ibm <- read.csv(paste0(workingDir_Results,"extn_ibm_MTN_1%.csv"))

## final population sizes
finalPop1 <- read.csv(paste0(workingDir_Results,"Nfinal_mtn0.99_IBM.csv"))
finalPop2 <- read.csv(paste0(workingDir_Results,"Nfinal_mtn0.85_IBM.csv"))
finalPop3 <- read.csv(paste0(workingDir_Results,"Nfinal_mtn0.65_IBM.csv"))
finalPop4 <- read.csv(paste0(workingDir_Results,"Nfinal_wlg0.42_IBM.csv"))
finalPop5 <- read.csv(paste0(workingDir_Results,"Nfinal_mtn3%_LM.csv"))
finalPop6 <- read.csv(paste0(workingDir_Results,"Nfinal_mtn2%_LM.csv"))
finalPop7 <- read.csv(paste0(workingDir_Results,"Nfinal_mtn1%_LM.csv"))
finalPop8 <- read.csv(paste0(workingDir_Results,"Nfinal_wlg_LM.csv"))

########################################
######## DEMOGRAPHIC PYRAMIDS ##########
########################################

## Demographic pyramid for MTN (same for 3%, 2%, and 1%)
n <- rep(1, nrow(dat))
n[1] <- 1
for (i in 2:length(n)){
  n[i] <- prod(1-dat[1:(i-1),2])
} ## Make sure sum equals 1 to generate pyramid
n_mtn <- n/(sum(n))
## for the cumulative survival curve:
n_mtnCS <- n/n[1]

## Demographic pyramid for WLG
n <- rep(1, nrow(dat))
n[1] <- 1
for (i in 2:length(n)){
  n[i] <- prod(1-dat[1:(i-1),4])
} ## Make sure sum equals 1 to generate pyramid
n_wlg <- n/(sum(n))
## for the cumulative survival curve:
n_wlgCS <- n/n[1]

## Select file name:
file_name <- "Fig2a_demographic_pyramid.pdf"

pdf(file_name, width=5,height=8)
## Dem Pyramid for MTN with 3%, 2% and 1% growth rate
barplot(n_mtn, horiz=T, names.arg=paste0(0:(length(n_mtn)-1), "-", 1:length(n_mtn)), las=1, xlab="Relative frequency", col=colMTN3, cex.axis = 1, cex.names = 0.7, ylab="Age", cex.lab=1, font.lab=2, xlim=c(0,0.07))

## Dem Pyramid for WLG
barplot(n_wlg, horiz=T, names.arg=paste0(0:(length(n_wlg)-1), "-", 1:length(n_wlg)), las=1, xlab="Relative frequency", col=coral, cex.axis = 1, cex.names = 0.7, ylab="Age", cex.lab=1, font.lab=2, add=TRUE)
dev.off()

########################################
###### CUMULATIVE SURVIVAL PLOTS #######
########################################

## Select file name:
file_name <- "Fig2b_survival_plots.pdf"

pdf(file_name, width=8,height=5)
par(mfrow=c(1,2), oma=c(4,1,1,1), mar=c(5,4,2,1))
plot(dat[,c(1,2)], type="o", bty="l", xlab="Age", ylab="Annual mortality", las=1, bg=colMTN3, cex.lab=0.8, cex.axis=0.8, ylim=c(0,1), cex=0.8, cex.lab=0.8, font.lab=2, pch=24)
lines(dat[,c(1,4)], type="o", bty="l", xlab="Age", ylab="Annual mortality", las=1, bg=colWLG, cex.lab=0.8, cex.axis=0.8, ylim=c(0,1), cex=0.8, cex.lab=0.8, font.lab=2, pch=21)
plot(dat[,1],n_mtnCS, bty="l", type="o", xlab="Age", ylab="Cumulative survival", las=1, bg=colMTN3, cex.lab=0.8, cex.axis=0.8, ylim=c(0,1), cex=0.8, cex.lab=0.8, font.lab=2, pch=24)
lines(dat[,1],n_wlgCS, bty="l", type="o", xlab="Age", ylab="Cumulative survival", las=1, bg=colWLG, cex.lab=0.8, cex.axis=0.8, ylim=c(0,1), cex=0.8, cex.lab=0.8, font.lab=2, pch=21)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("Mountain Gorillas", "Western Lowland Gorillas"),
       pt.bg=c(colMTN3, colWLG), lty=c(1,1), cex=0.8, text.font=2, pch=c(24,21), xpd = TRUE, horiz = FALSE, inset = c(0, 0.05), bty = "y")
dev.off()

# xpd = TRUE tells R that it is OK to plot outside the region 
# inset = c(x,y) tells R how to move the legend relative to the 'bottom' location

#########################################
######## EXTINCTION RISK PLOTS  #########
#########################################

## Select file name:
file_name <- "Fig4_extn_risk_LM.pdf"

pdf(file_name, width=6,height=6)
plot.new()
par(mfrow=c(1,1), oma=c(4,1,1,1), mar=c(5,4,2,1))
plot.window(xlim=c(1,9), ylim=c(0,100))
axis(1, 1:9, LETTERS[1:9])
axis(2)
axis(2, font.lab=2, at=seq(0, 100, by=10), labels=seq(0, 100, by=10))
title(xlab="Reintroduction Scenario", ylab="Extinction Risk (LM)", font.lab=2)
lines(wlg_lm$scenario, wlg_lm$prob_Extn, bg="#FF7F50", type="b", pch=21)
lines(mtn_1per_lm$scenario, mtn_1per_lm$prob_Extn, bg="#D58363", type="b", pch=22)
lines(mtn_2per_lm$scenario, mtn_2per_lm$prob_Extn, bg="#AC8777", type="b", pch=23)
lines(mtn_3per_lm$scenario, mtn_3per_lm$prob_Extn, bg="#838B8B", type="b", pch=24)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("Western Gorillas", "Mountain Gorillas - 1%", "Mountain Gorillas - 2%", "Mountain Gorillas - 3%"),
       pt.bg=c("#FF7F50", "#D58363","#AC8777","#838B8B"), lty=c(1,1,1,1), text.font=2, pch=c(21,22,23,24), cex=0.6, xpd = TRUE, horiz = FALSE, inset = c(0, 0.01), bty = "y")
dev.off()

## Select file name:
file_name <- "Fig4_extn_risk_IBM.pdf"

pdf(file_name, width=6,height=6)
plot.new()
par(mfrow=c(1,1), oma=c(4,1,1,1), mar=c(5,4,2,1))
plot.window(xlim=c(1,9), ylim=c(0,100))
axis(1, 1:9, LETTERS[1:9])
axis(2)
axis(2, font.lab=2, at=seq(0, 100, by=10), labels=seq(0, 100, by=10))
title(xlab="Reintroduction Scenario", ylab="Extinction Risk (IBM)", font.lab=2)
lines(wlg_ibm$scenario, wlg_ibm$prob_Extn, bg="#FF7F50", type="b", pch=21, lty=2)
lines(mtn_1per_ibm$scenario, mtn_1per_ibm$prob_Extn, bg="#D58363", type="b", pch=22, lty=2)
lines(mtn_2per_ibm$scenario, mtn_2per_ibm$prob_Extn, bg="#AC8777", type="b", pch=23, lty=2)
lines(mtn_3per_ibm$scenario, mtn_3per_ibm$prob_Extn, bg="#838B8B", type="b", pch=24, lty=2)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("Western Gorillas", "Mountain Gorillas - 1%", "Mountain Gorillas - 2%", "Mountain Gorillas - 3%"),
       pt.bg=c("#FF7F50", "#D58363","#AC8777","#838B8B"), lty=c(2,2,2,2), text.font=2, pch=c(21,22,23,24), cex=0.6, xpd = TRUE, horiz = FALSE, inset = c(0, 0.01), bty = "y")
dev.off()

########################################
#### LESLIE MATRIX PROJECTION PLOTS ####
########################################

## Select file name:
file_name <- "FigS1_LM_Projections_MTN3.pdf"
##file_name <- "FigS1_LM_Projections_MTN2.pdf"
##file_name <- "FigS1_LM_Projections_MTN1.pdf"
##file_name <- "FigS1_LM_Projections_WLG.pdf"

## Select exticntion risk files:
##probExt_lm <- read.csv("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_extn_results/extn_lm_MTN_3%.csv")
##probExt_lm <- read.csv("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_extn_results/extn_lm_MTN_2%.csv")
##probExt_lm <- read.csv("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_extn_results/extn_lm_MTN_1%.csv")
##probExt_lm <- read.csv("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_extn_results/extn_lm_WLG.csv")

## Select the correct folder for either WLG or MTN data
workingDir_LM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_LM_50year_mtn_3%/")
##workingDir_LM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_LM_50year_mtn_2%/")
##workingDir_LM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_LM_50year_mtn_1%/")
##workingDir_LM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_LM_50year_wlg/")

setwd(workingDir_LM)
allScenarioFiles <- list.files(pattern="*.csv")

for (i in 1:length(allScenarioFiles)){
  assign(allScenarioFiles[i], 
         read.csv(paste(workingDir_LM, allScenarioFiles[i], sep=''), header=TRUE)
  )}

stochObjects <- c("LM_Scenario1.csv","LM_Scenario2.csv","LM_Scenario3.csv","LM_Scenario4.csv","LM_Scenario5.csv","LM_Scenario6.csv","LM_Scenario7.csv","LM_Scenario8.csv","LM_Scenario9.csv")
detObjects <- c("LM_Det_Scenario1.csv","LM_Det_Scenario2.csv","LM_Det_Scenario3.csv","LM_Det_Scenario4.csv", "LM_Det_Scenario5.csv","LM_Det_Scenario6.csv","LM_Det_Scenario7.csv","LM_Det_Scenario8.csv","LM_Det_Scenario9.csv")

setwd(workingDir)
pdf(file_name, width=10,height=10)
par(mfrow=c(3,3), oma=c(0,0,0,0), mar=c(5,4,2,1), las=1, bty="l")
maxY <- 120 ## max y-axis value
time <- 0:50 ## time interval for the plots
for(j in 1:length(stochObjects)){
  probExt <- probExt_lm[j,3]
  scenario <- as.character(probExt_lm[j,1])
  tempX <- get(stochObjects[j])
  N_projected_detX <- unname(unlist(get(detObjects[j])))
  plot(N_projected_detX~time, type="l", col=1, xlab="Years post-introduction", ylab="Population Size", ylim=c(0,maxY),lty=2, cex.lab=1, cex.axis=1, font.lab=2) ## plot of deterministic projection
  for(i in 1:ncol(tempX)){
    lines(time, tempX[,i], col=grey(.8, alpha=.05), lwd=3)
  } ## plots projections from stochastic LM simulations
  lines(apply(tempX, 1, mean)~time, type="l", col=2, lwd=4) ## plot mean projection from stochastic LM simulations
  lines(N_projected_detX~time, type="l", col=1, lwd=2, lty=2) ## replot deterministic projection
  title(main=paste0("Scenario ",scenario), sub=paste0("Extinction Risk = ",probExt, "%"), cex.main=1, cex.sub=1, col.sub=1, font.sub=3)
  qtiles <- apply(tempX, 1, function(v) quantile(v, probs=c(0.05, 0.95))) ## plot 95% confidence intervals for simulations
  lines((0:(nrow(tempX)-1)), qtiles[1,], col=1, lty=2)
  lines((0:(nrow(tempX)-1)), qtiles[2,], col=1, lty=2)
  lines(x=c(-5:50), y=rep(50, 56), col="navyblue", lwd=2, lty=1) ## add a line for the 50 individual mark
}
dev.off()

########################################
######### IBM PROJECTION PLOTS #########
########################################

## Select file name:
file_name <- "FigS2_IBM_Projections_MTN3.pdf"
##file_name <- "FigS2_IBM_Projections_MTN2.pdf"
##file_name <- "FigS2_IBM_Projections_MTN1.pdf"
##file_name <- "FigS2_IBM_Projections_WLG.pdf"


## Select exticntion risk files:
probExt_ibm <- read.csv("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_extn_results/extn_ibm_MTN_3%.csv")
##probExt_ibm <- read.csv("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_extn_results/extn_ibm_MTN_2%.csv")
##probExt_ibm <- read.csv("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_extn_results/extn_ibm_MTN_1%.csv")
##probExt_ibm <- read.csv("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_extn_results/extn_ibm_WLG.csv")

## Select the correct folder for either WLG or MTN data
workingDir_IBM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_IBM_50year_mtn_0.99/")
##workingDir_IBM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_IBM_50year_mtn_0.85/")
##workingDir_IBM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_IBM_50year_mtn_0.65/")
##workingDir_IBM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_IBM_50year_wlg_0.42/")

setwd(workingDir_IBM)
allScenarioFiles <- list.files(pattern="*.csv")

for (i in 1:length(allScenarioFiles)){
  assign(allScenarioFiles[i], 
         read.csv(paste(workingDir_IBM, allScenarioFiles[i], sep=''), header=TRUE)
  )}

stochObjects <- c("IBM_Scenario1.csv","IBM_Scenario2.csv","IBM_Scenario3.csv","IBM_Scenario4.csv","IBM_Scenario5.csv","IBM_Scenario6.csv","IBM_Scenario7.csv","IBM_Scenario8.csv","IBM_Scenario9.csv")

timeunit <- 1/12 ## time interval for the plots

setwd(workingDir)
pdf(file_name, width=10,height=10) 
par(mfrow=c(3,3), oma=c(0,0,0,0), mar=c(5,4,2,1), las=1, bty="l")
maxY <- 120 ## max y-axis value
for(j in 1:length(stochObjects)){
  probExt <- probExt_ibm[j,3]
  scenario <- as.character(probExt_ibm[j,1])
  resX <- get(stochObjects[j])
  
  plot((0:(nrow(resX)-1))*timeunit, apply(resX, 1, mean), type="l", col=2, lwd=2, xlab="Years post-introduction", ylab="Population size", font.lab=2, bty="l", ylim=c(0,120))
  for(i in 1:ncol(resX)){
    lines((0:(nrow(resX)-1))*timeunit, resX[,i], col=grey(.8, alpha=.05), lwd=3)
  }
  
  ## add mean trend
  lines((0:(nrow(resX)-1))*timeunit, apply(resX, 1, mean), type="l", col=2, lwd=3)
  
  ## add 95% upper/lower limits
  qtiles <- apply(resX, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
  lines((0:(nrow(resX)-1))*timeunit, qtiles[1,], lty=2)
  lines((0:(nrow(resX)-1))*timeunit, qtiles[2,], lty=2)
  
  title(main=paste0("Scenario ",scenario), sub=paste0("Extinction Risk = ",probExt, "%"), cex.main=1, cex.sub=1, col.sub=1,  font.sub=3)
  
  ## add a line for the 50 individual mark
  lines(x=c(-5:50), y=rep(50, 56), col="navy", lwd=2, lty=1)
}
dev.off()

########################################
######## FINAL POPULATION SIZES ########
########################################
setwd(workingDir)

## Select file name:
file_name <- "Fig5_FinalPopSize_IBM.pdf"

pdf(file_name, width=8,height=8)
par(mfrow=c(2,2))
boxplot(finalPop1$A, finalPop1$B, finalPop1$C, finalPop1$D, finalPop1$E, finalPop1$F, finalPop1$G, finalPop1$H, finalPop1$I, pch=24, col=colMTN3, bg=colMTN3, ylim=c(0,300), varwidth=FALSE)
lines(x=c(-5:51), y=rep(50, 57), col="red", lwd=2, lty=2)
par(new = TRUE, mar=c(5.1,4.1,4.1,5.1))
axis(2, at=seq(0, 150, by=50), labels=seq(0, 150, by=50), ylab="Population size after 50 years", font.lab=2)
axis(1, LETTERS[1:9], at=1:9, labels=LETTERS[1:9])
title(xlab="Reintroduction Scenario", ylab="Population size after 50 years", font.lab=2)

boxplot(finalPop2$A, finalPop2$B, finalPop2$C, finalPop2$D, finalPop2$E, finalPop2$F, finalPop2$G, finalPop2$H, finalPop2$I, pch=23, col=colMTN2, bg=colMTN2, ylim=c(0,300), varwidth=FALSE)
lines(x=c(-5:51), y=rep(50, 57), col="red", lwd=2, lty=2)
par(new = TRUE, mar=c(5.1,4.1,4.1,5.1))
axis(2, at=seq(0, 150, by=50), labels=seq(0, 150, by=50), ylab="Population size after 50 years", font.lab=2)
axis(1, LETTERS[1:9], at=1:9, labels=LETTERS[1:9])
title(xlab="Reintroduction Scenario", ylab="Population size after 50 years", font.lab=2)

boxplot(finalPop3$A, finalPop3$B, finalPop3$C, finalPop3$D, finalPop3$E, finalPop3$F, finalPop3$G, finalPop3$H, finalPop3$I, pch=22, col=colMTN2, bg=colMTN2, ylim=c(0,300), varwidth=FALSE)
lines(x=c(-5:51), y=rep(50, 57), col="red", lwd=2, lty=2)
par(new = TRUE, mar=c(5.1,4.1,4.1,5.1))
axis(2, at=seq(0, 150, by=50), labels=seq(0, 150, by=50), ylab="Population size after 50 years", font.lab=2)
axis(1, LETTERS[1:9], at=1:9, labels=LETTERS[1:9])
title(xlab="Reintroduction Scenario", ylab="Population size after 50 years", font.lab=2)

boxplot(finalPop4$A, finalPop4$B, finalPop4$C, finalPop4$D, finalPop4$E, finalPop4$F, finalPop4$G, finalPop4$H, finalPop4$I, pch=21, col=colMTN2, bg=colMTN2, ylim=c(0,300), varwidth=FALSE)
lines(x=c(-5:51), y=rep(50, 57), col="red", lwd=2, lty=2)
par(new = TRUE, mar=c(5.1,4.1,4.1,5.1))
axis(2, at=seq(0, 150, by=50), labels=seq(0, 150, by=50), ylab="Population size after 50 years", font.lab=2)
axis(1, LETTERS[1:9], at=1:9, labels=LETTERS[1:9])
title(xlab="Reintroduction Scenario", ylab="Population size after 50 years", font.lab=2)
dev.off()

## Select file name:
file_name <- "Fig5_FinalPopSize_LM.pdf"

pdf(file_name, width=8,height=8)
par(mfrow=c(2,2))
boxplot(finalPop5$A, finalPop5$B, finalPop5$C, finalPop5$D, finalPop5$E, finalPop5$F, finalPop5$G, finalPop5$H, finalPop5$I, pch=24, col=colMTN3, bg=colMTN3, ylim=c(0,300), varwidth=FALSE)
lines(x=c(-5:51), y=rep(50, 57), col="red", lwd=2, lty=2)
par(new = TRUE, mar=c(5.1,4.1,4.1,5.1))
axis(2, at=seq(0, 150, by=50), labels=seq(0, 150, by=50), ylab="Population size after 50 years", font.lab=2)
axis(1, LETTERS[1:9], at=1:9, labels=LETTERS[1:9])
title(xlab="Reintroduction Scenario", ylab="Population size after 50 years", font.lab=2)

boxplot(finalPop6$A, finalPop6$B, finalPop6$C, finalPop6$D, finalPop6$E, finalPop6$F, finalPop6$G, finalPop6$H, finalPop6$I, pch=23, col=colMTN2, bg=colMTN2, ylim=c(0,300), varwidth=FALSE)
lines(x=c(-5:51), y=rep(50, 57), col="red", lwd=2, lty=2)
par(new = TRUE, mar=c(5.1,4.1,4.1,5.1))
axis(2, at=seq(0, 150, by=50), labels=seq(0, 150, by=50), ylab="Population size after 50 years", font.lab=2)
axis(1, LETTERS[1:9], at=1:9, labels=LETTERS[1:9])
title(xlab="Reintroduction Scenario", ylab="Population size after 50 years", font.lab=2)

boxplot(finalPop7$A, finalPop7$B, finalPop7$C, finalPop7$D, finalPop7$E, finalPop7$F, finalPop7$G, finalPop7$H, finalPop7$I, pch=22, col=colMTN2, bg=colMTN2, ylim=c(0,300), varwidth=FALSE)
lines(x=c(-5:51), y=rep(50, 57), col="red", lwd=2, lty=2)
par(new = TRUE, mar=c(5.1,4.1,4.1,5.1))
axis(2, at=seq(0, 150, by=50), labels=seq(0, 150, by=50), ylab="Population size after 50 years", font.lab=2)
axis(1, LETTERS[1:9], at=1:9, labels=LETTERS[1:9])
title(xlab="Reintroduction Scenario", ylab="Population size after 50 years", font.lab=2)

boxplot(finalPop8$A, finalPop8$B, finalPop8$C, finalPop8$D, finalPop8$E, finalPop8$F, finalPop8$G, finalPop8$H, finalPop8$I, pch=21, col=colMTN2, bg=colMTN2, ylim=c(0,300), varwidth=FALSE)
lines(x=c(-5:51), y=rep(50, 57), col="red", lwd=2, lty=2)
par(new = TRUE, mar=c(5.1,4.1,4.1,5.1))
axis(2, at=seq(0, 150, by=50), labels=seq(0, 150, by=50), ylab="Population size after 50 years", font.lab=2)
axis(1, LETTERS[1:9], at=1:9, labels=LETTERS[1:9])
title(xlab="Reintroduction Scenario", ylab="Population size after 50 years", font.lab=2)
dev.off()

#####################################################################################
################ Historical population trends for Tshiaberimu gorillas ##############
#####################################################################################

## Select file name:
file_name <- "Fig1_HistoricalTshiaberimuPop.pdf"

## Take a look at the historical population trajectories using Tshiaberimu census data
year <- c(1959,1986,1995,1996,2003,2004,2006,2007,2008,2009,2011,2012,2013,2016,2017) ## census years
numyears <- 1959:2017 ## census time period
N <- c(35,20,17,16,20,20,21,22,18,16,6,6,7,6,6) ## census data

## let's look at the rate of change in this population
## lambda is the finite rate of increase of a population over one time step. r is the intrinsinc rate of growth. negative r values indicate a population in decline. lambda < 1 indicates a decline. the relationship between lambda and r : lambda = Nt+1  / Nt, r = ln(lambda), lambda = e^r
logLambda <- (1/58)*log(6/35) ## 58 years for the census time period, loglambda = 1/timeperiod*log(Ntfinal)/Nt0
lambda <- exp(logLambda)
popEst <- 35*(exp(logLambda))^(0:58) ## this is the expected rate of change in the population given Ntfinal and Nt0
popEst ## these are the predicted population estimates given the calculated lambda value

## let's fit these parameter estimates to a linear model to calculate the r and lambda values to get a more accurate estimate of these parameters:
modelGeom <- lm(log(N)~year) ## should be linear on a log scale
r_lm <- modelGeom$coef[2] ## take the slope of the line from this linear model for the intrinsic rate of growth r=-0.0289175
lambda_lm <- exp(modelGeom$coef[2]) ## lambda=0.9714966
popEst_lm <- 35*(exp(r_lm))^(0:58)

pdf(file_name, width=5,height=5)
## plot the actual population sizes from census data and the expected population size:
plot(year, N, xlab="Census Year", 
     ylab="Estimated population size, N", 
     pch=19, type="o",
     ylim=c(0,50), 
     xlim=c(numyears[1],numyears[length(numyears)]), 
     las=1, cex.main=0.8, cex.lab=0.8, cex.axis=0.8, font.lab=2)
lines(numyears, popEst_lm, col=2, lty=2, lwd=2)
dev.off()
