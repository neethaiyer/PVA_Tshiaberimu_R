## Set the working directory:
workingDir <- "~/Box Sync/PVA_Paper/PVA_Tshiaberimu_R/"
setwd(workingDir)
workingDir_Results <- "~/Box Sync/PVA_Paper/PVA_Tshiaberimu_R/pva_extn_results/"
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

## Dem Pyramid for MTN with 3%, 2% and 1% growth rate
barplot(n_mtn, horiz=T, names.arg=paste0(0:(length(n_mtn)-1), "-", 1:length(n_mtn)), las=1, xlab="relative frequency", col=colMTN3, cex.axis = 1, cex.names = 0.7, ylab="Age", cex.lab=1, font.lab=2, xlim=c(0,0.07))

## Dem Pyramid for WLG
barplot(n_wlg, horiz=T, names.arg=paste0(0:(length(n_wlg)-1), "-", 1:length(n_wlg)), las=1, xlab="relative frequency", col=coral, cex.axis = 1, cex.names = 0.7, ylab="Age", cex.lab=1, font.lab=2, add=TRUE)

########################################
###### CUMULATIVE SURVIVAL PLOTS #######
########################################

par(mfrow=c(1,2), oma=c(4,1,1,1), mar=c(5,4,2,1))
plot(dat[,c(1,2)], type="o", bty="l", xlab="Age", ylab="Annual mortality", las=1, bg=colMTN3, cex.lab=0.8, cex.axis=0.8, ylim=c(0,1), cex=0.8, cex.lab=0.8, font.lab=2, pch=24)
lines(dat[,c(1,4)], type="o", bty="l", xlab="Age", ylab="Annual mortality", las=1, bg=colWLG, cex.lab=0.8, cex.axis=0.8, ylim=c(0,1), cex=0.8, cex.lab=0.8, font.lab=2, pch=21)
plot(dat[,1],n_mtnCS, bty="l", type="o", xlab="Age", ylab="Cumulative survival", las=1, bg=colMTN3, cex.lab=0.8, cex.axis=0.8, ylim=c(0,1), cex=0.8, cex.lab=0.8, font.lab=2, pch=24)
lines(dat[,1],n_wlgCS, bty="l", type="o", xlab="Age", ylab="Cumulative survival", las=1, bg=colWLG, cex.lab=0.8, cex.axis=0.8, ylim=c(0,1), cex=0.8, cex.lab=0.8, font.lab=2, pch=21)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("Mountain Gorillas", "Western Lowland Gorillas"),
       pt.bg=c(colMTN3, colWLG), lty=c(1,1), cex=0.8, text.font=2, pch=c(24,21), xpd = TRUE, horiz = FALSE, inset = c(0, 0.05), bty = "y")

# xpd = TRUE tells R that it is OK to plot outside the region 
# inset = c(x,y) tells R how to move the legend relative to the 'bottom' location

#########################################
##### EXTINCTION PROBABILITY PLOTS  #####
#########################################

plot.new()
par(mfrow=c(1,1), oma=c(4,1,1,1), mar=c(5,4,2,1))
plot.window(xlim=c(1,9), ylim=c(0,100))
axis(1, 1:9, LETTERS[1:9])
axis(2)
axis(2, font.lab=2, at=seq(0, 100, by=10), labels=seq(0, 100, by=10))
title(xlab="Reintroduction Scenario", ylab="Likelihood of Extinction (LM)", font.lab=2)
lines(wlg_lm$scenario, wlg_lm$prob_Extn, bg="#FF7F50", type="b", pch=21)
lines(mtn_1per_lm$scenario, mtn_1per_lm$prob_Extn, bg="#D58363", type="b", pch=22)
lines(mtn_2per_lm$scenario, mtn_2per_lm$prob_Extn, bg="#AC8777", type="b", pch=23)
lines(mtn_3per_lm$scenario, mtn_3per_lm$prob_Extn, bg="#838B8B", type="b", pch=24)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("Western Gorillas", "Mountain Gorillas - 1%", "Mountain Gorillas - 2%", "Mountain Gorillas - 3%"),
       pt.bg=c("#FF7F50", "#D58363","#AC8777","#838B8B"), lty=c(1,1,1,1), text.font=2, pch=c(21,22,23,24), cex=0.6, xpd = TRUE, horiz = FALSE, inset = c(0, 0.01), bty = "y")

plot.new()
par(mfrow=c(1,1), oma=c(4,1,1,1), mar=c(5,4,2,1))
plot.window(xlim=c(1,9), ylim=c(0,100))
axis(1, 1:9, LETTERS[1:9])
axis(2)
axis(2, font.lab=2, at=seq(0, 100, by=10), labels=seq(0, 100, by=10))
title(xlab="Reintroduction Scenario", ylab="Likelihood of Extinction (IBM)", font.lab=2)
lines(wlg_ibm$scenario, wlg_ibm$prob_Extn, bg="#FF7F50", type="b", pch=21, lty=2)
lines(mtn_1per_ibm$scenario, mtn_1per_ibm$prob_Extn, bg="#D58363", type="b", pch=22, lty=2)
lines(mtn_2per_ibm$scenario, mtn_2per_ibm$prob_Extn, bg="#AC8777", type="b", pch=23, lty=2)
lines(mtn_3per_ibm$scenario, mtn_3per_ibm$prob_Extn, bg="#838B8B", type="b", pch=24, lty=2)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("Western Gorillas", "Mountain Gorillas - 1%", "Mountain Gorillas - 2%", "Mountain Gorillas - 3%"),
       pt.bg=c("#FF7F50", "#D58363","#AC8777","#838B8B"), lty=c(2,2,2,2), text.font=2, pch=c(21,22,23,24), cex=0.6, xpd = TRUE, horiz = FALSE, inset = c(0, 0.01), bty = "y")

########################################
#### LESLIE MATRIX PROJECTION PLOTS ####
########################################

simObjects <- c("tempA", "tempB", "tempC", "tempD", "tempE", "tempF", "tempG", "tempH", "tempI")
detObjects <- c("N_projected_detA", "N_projected_detB", "N_projected_detC", "N_projected_detD", "N_projected_detE", "N_projected_detF", "N_projected_detG", "N_projected_detH", "N_projected_detI")
probExt_lm <- prob_50years_mtn_lm

par(mfrow=c(3,3), oma=c(0,0,0,0), mar=c(5,4,2,1), las=1, bty="l")
maxY <- 120 ## max y-axis value
time <- 0:50 ## time interval for the plots

for(j in 1:9){
  ##probExt <- prob_50years_wlg_lm[j,6] ##uncomment this line for western gorilla demographic parameter simulations
  probExt <- prob_50years_mtn_lm[j,6] ##uncomment this line for mountain gorilla demographic parameter simulations
  tempX <- get(simObjects[j])
  N_projected_detX <- get(detObjects[j])
  plot(N_projected_detX~time, type="l", col=1, xlab="Years post-introduction", ylab="Population Size", ylim=c(0,maxY),lty=2, cex.lab=1, cex.axis=1, font.lab=2) ## plot of deterministic projection
  for(i in 1:ncol(tempX)){
    lines(time, tempX[,i], col=grey(.8, alpha=.05), lwd=3)
  } ## plots projections from stochastic LM simulations
  lines(apply(tempX, 1, mean)~time, type="l", col=2, lwd=4) ## plot mean projection from stochastic LM simulations
  lines(N_projected_detX~time, type="l", col=1, lwd=2, lty=2) ## replot deterministic projection
  title(main="A: No introduction, male juvenile", sub=paste0("Probability of extinction = ",probExt, "%"), cex.main=1, cex.sub=1, col.sub=1, font.sub=3)
  qtiles <- apply(tempX, 1, function(v) quantile(v, probs=c(0.05, 0.95))) ## plot 95% confidence intervals for simulations
  lines((0:(nrow(tempX)-1)), qtiles[1,], col=1, lty=2)
  lines((0:(nrow(tempX)-1)), qtiles[2,], col=1, lty=2)
  lines(x=c(-5:50), y=rep(50, 56), col="navyblue", lwd=2, lty=1) ## add a line for the 50 individual mark
}
