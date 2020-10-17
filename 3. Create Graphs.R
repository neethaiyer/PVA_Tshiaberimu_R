## Set the working directory:
workingDir <- "/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/"
setwd(workingDir)
workingDir_Output <- "/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/PVA_Output/"
workingDir_Input <- "/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/PVA_Input/"
workingDir_Figures <- "/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/PVA_Figures/"
##workingDir <- "~/Documents/git repositories/PVA_Tshiaberimu_R/"

## Let's pick our colors for MTN and WLG and include a transparency factor; note that these are only used for the bargraph below:
coral <- rgb(250, 128, 114, alpha=150, maxColorValue = 255) ## for WLG
azure4 <- rgb(131, 139, 139, alpha=150, maxColorValue = 255) ## for MTN
colfun <- colorRampPalette(c("azure4", "coral"))
colfun(4)
colMTN3_alpha <- rgb(37, 37, 37, alpha=150, maxColorValue = 255)
colMTN2_alpha <- rgb(99, 99, 99, alpha=150, maxColorValue = 255)
colMTN1_alpha <- rgb(150, 150, 150, alpha=150, maxColorValue = 255)
colWLG_alpha <- rgb(189, 189, 189, alpha=150, maxColorValue = 255)

## A. Color for MTN with 3% growth rate
##colMTN3 <- "#838B8B"
colMTN3 <- "#252525"
## B. Color for MTN with 2% growth rate
##colMTN2 <- "#AC8777"
colMTN2 <- "#636363"
## C. Color for MTN with 1% growth rate
##colMTN1 <- "#D58363"
colMTN1 <- "#969696"
## D. Color for WLG 
##colWLG <- "#FF7F50"
colWLG <- "#bdbdbd"
## Color for N=50 mark
colN50 <- "#de2d26"

##plot(rep(1,5), col=c(colMTN3,colMTN2,colMTN1,colWLG, colN50), pch=19 ,cex=3)

## Read all csv files. Your life history tables should have at least 3 columns: age, mortality rate, and fertility rate:
dat <- read.csv(paste0(workingDir_Input, "Gorilla_LifeTables.csv"))
dat$fertilityrate_MTN1 <- dat[,3]*.789 
## fertility rates multiplied by factor less than 1 to get eigen values of 1.01 which corresponds to a 1% growth rate
dat$fertilityrate_MTN2 <- dat[,3]*.643 
## fertility rates multiplied by factor less than 1 to get eigen values of 1.02 which corresponds to a 2% growth rate

## probability of extinctions based on LM projections
mtn_3per_lm <- read.csv(paste0(workingDir_Output,"Results/Results_LM_mtn_3%.csv"))
mtn_2per_lm <- read.csv(paste0(workingDir_Output,"Results/Results_LM_mtn_2%.csv"))
mtn_1per_lm <- read.csv(paste0(workingDir_Output,"Results/Results_LM_mtn_1%.csv"))
##wlg_lm <- read.csv(paste0(workingDir_Output,"Results/Results_LM_wlg.csv"))
mtn_3per_ibm <- read.csv(paste0(workingDir_Output,"Results/Results_IBM_mtn_3%.csv"))
mtn_2per_ibm <- read.csv(paste0(workingDir_Output,"Results/Results_IBM_mtn_2%.csv"))
mtn_1per_ibm <- read.csv(paste0(workingDir_Output,"Results/Results_IBM_mtn_1%.csv"))
##wlg_ibm <- read.csv(paste0(workingDir_Output,"Results/Results_IBM_wlg.csv"))

## final population sizes
finalPop1 <- read.csv(paste0(workingDir_Output,"Results/Results_LM_Nfinal_mtn_3%.csv"))
finalPop2 <- read.csv(paste0(workingDir_Output,"Results/Results_LM_Nfinal_mtn_2%.csv"))
finalPop3 <- read.csv(paste0(workingDir_Output,"Results/Results_LM_Nfinal_mtn_1%.csv"))
##finalPop4 <- read.csv(paste0(workingDir_Output,"Results/Results_LM_Nfinal_wlg.csv"))
finalPop5 <- read.csv(paste0(workingDir_Output,"Results/Results_IBM_Nfinal_mtn_3%.csv"))
finalPop6 <- read.csv(paste0(workingDir_Output,"Results/Results_IBM_Nfinal_mtn_2%.csv"))
finalPop7 <- read.csv(paste0(workingDir_Output,"Results/Results_IBM_Nfinal_mtn_1%.csv"))
##finalPop8 <- read.csv(paste0(workingDir_Output,"Results/Results_IBM_Nfinal_wlg.csv"))


#####################################################################################
########## FIGURE 1: Historical population trends for Tshiaberimu gorillas ##########
#####################################################################################

## Select file name:
file_name <- "Fig1_HistoricalTshiaberimuPop.pdf"

## Take a look at the historical population trajectories using Tshiaberimu census data
year <- c(1959,1986,1995,1996,2003,2004,2006,2007,2008,2009,2011,2012,2013,2015,2017,2020) ## census years
censusPeriod <- 1959:2020 ## vector for census period
totalYears <- length(censusPeriod)-1
N <- c(35,20,17,16,20,20,21,22,18,16,6,6,7,6,6,6) ## census data

## let's look at the rate of change in this population
## lambda is the finite rate of increase of a population over one time step. r is the intrinsinc rate of growth. negative r values indicate a population in decline. lambda < 1 indicates a decline. the relationship between lambda and r : lambda = Nt+1  / Nt, r = ln(lambda), lambda = e^r
logLambda <- (1/totalYears)*log(N[length(N)]/N[1]) ## loglambda = 1/timeperiod*log(Ntfinal)/Nt0
lambda <- exp(logLambda)
popEst <- N[1]*(exp(logLambda))^(0:totalYears) ## this is the expected rate of change in the population given Ntfinal and Nt0
popEst ## these are the predicted population estimates given the calculated lambda value

## let's fit these parameter estimates to a linear model to calculate the r and lambda values to get a more accurate estimate of these parameters:
modelGeom <- lm(log(N)~year) ## should be linear on a log scale
r_lm <- modelGeom$coef[2] ## take the slope of the line from this linear model for the intrinsic rate of growth r=-0.0297
lambda_lm <- exp(modelGeom$coef[2]) ## lambda=0.97
popEst_lm <- N[1]*(exp(r_lm))^(0:totalYears)

setwd(workingDir_Figures)
pdf(file_name, width=5,height=5)
## plot the actual population sizes from census data and the expected population size:
plot(year, N, 
     xlab="Census Year", ylab="Population Size", 
     pch=16, type="o",
     ylim=c(0,40), 
     cex.lab=0.8, cex.axis=0.8, font.lab=2, axes=FALSE)
box(bty="l")
axis(2, font.lab=2, at=seq(0, 50, by=10), labels=seq(0, 50, by=10), las=2)
axis(1, font.lab=2, at=seq(censusPeriod[1]+1, censusPeriod[length(censusPeriod)]+1, by=10), labels=seq(censusPeriod[1]+1, censusPeriod[length(censusPeriod)]+1, by=10), las=1)
lines(censusPeriod, popEst_lm, col=colN50, lty=2, lwd=2)
dev.off()

#####################################################################################
######################## FIGURE 2: CUMULATIVE SURVIVAL PLOTS ########################
#####################################################################################

## Eigenvector for mountain gorillas (same for 3%, 2%, and 1%)
n <- rep(1, nrow(dat))
n[1] <- 1
for (i in 2:length(n)){
  n[i] <- prod(1-dat[1:(i-1),2])
} ## Make sure sum equals 1 to generate pyramid
n_mtn <- n/(sum(n))
## for the cumulative survival curve:
n_mtnCS <- n/n[1]

## Eigenvector for wester lowland gorillas
n <- rep(1, nrow(dat))
n[1] <- 1
for (i in 2:length(n)){
  n[i] <- prod(1-dat[1:(i-1),4])
} ## Make sure sum equals 1 to generate pyramid
n_wlg <- n/(sum(n))
## for the cumulative survival curve:
n_wlgCS <- n/n[1]

## Select file name:
file_name <- "Fig2_survival_plot.pdf"

setwd(workingDir_Figures)
pdf(file_name, width=5,height=5)
par(oma=c(4,1,1,1), mar=c(5,4,2,1))
plot(dat[,1],n_mtnCS, bty="l", type="o", xlab="Age (years)", ylab="Cumulative survival rate", las=1, bg=colMTN3, cex.lab=0.8, cex.axis=0.8, ylim=c(0,1), cex=0.7, cex.lab=0.8, font.lab=2, pch=24)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
dev.off()

#####################################################################################
######################### FIGURE 4: EXTINCTION RISK PLOTS  ##########################
#####################################################################################

## Select file name:
file_name <- "Fig4_extn_risk_LM_IBM.pdf"

setwd(workingDir_Figures)
pdf(file_name, width=8,height=4)
layout(matrix(c(1,2,3,4,4,4), ncol=3, byrow=TRUE), heights=c(0.8,0.2))
par(oma=c(1, 1, 1, 1), mar=c(4, 4, 2, 2))

plot(mtn_1per_lm$extn_Risk, bg=colMTN1, type="o", pch=21, xlab=NA, ylab=NA, xlim=c(1,9), ylim=c(0,50), xaxt="n", yaxt = "n", axes=FALSE)
lines(mtn_1per_ibm$scenario, mtn_1per_ibm$extn_Risk, bg=colMTN1, type="b", pch=21, lty=2)
box(bty="l")
axis(1, font.lab=2, at=seq(1, 9, by=1), labels=LETTERS[1:9], las=1)
axis(2, font.lab=2, at=seq(0, 50, by=10), labels=seq(0, 50, by=10), las=2)
title(main="Mountain gorilla, 1%", ylab="Extinction risk (%)", cex=0.8, cex.main=0.9, font.lab=2)

plot(mtn_2per_lm$extn_Risk, bg=colMTN2, type="o", pch=23, xlab=NA, ylab=NA, xlim=c(1,9), ylim=c(0,50), xaxt="n", yaxt = "n", axes=FALSE)
lines(mtn_2per_ibm$scenario, mtn_2per_ibm$extn_Risk, bg=colMTN2, type="b", pch=23, lty=2)
box(bty="l")
axis(1, font.lab=2, at=seq(1, 9, by=1), labels=LETTERS[1:9], las=1)
axis(2, font.lab=2, at=seq(0, 50, by=10), labels=seq(0, 50, by=10), las=2)
title(main="Mountain gorilla, 2%", cex.main=0.9, font.lab=2)

plot(mtn_3per_lm$extn_Risk, bg=colMTN3, type="o", pch=24, xlab=NA, ylab=NA, xlim=c(1,9), ylim=c(0,50), xaxt="n", yaxt = "n", axes=FALSE)
lines(mtn_3per_ibm$extn_Risk, bg=colMTN3, type="b", pch=24, lty=2)
box(bty="l")
axis(1, font.lab=2, at=seq(1, 9, by=1), labels=LETTERS[1:9], las=1)
axis(2, font.lab=2, at=seq(0, 50, by=10), labels=seq(0, 50, by=10), las=2)
title(main="Mountain gorilla, 3.2%", cex.main=0.9, font.lab=2)

mtext("Reintroduction Scenario", side=1, line=-6, outer=TRUE, cex=0.7, font=2)

par(mai=c(0,0,0,0))
plot.new()
legend("bottom", legend=c("Leslie Matrix", "Individual-based Model"), lty=c(1,2), text.font=2, xpd = TRUE, horiz = FALSE, inset = c(0, 0.04), bty = "y")
dev.off()

# xpd = TRUE: legend will go outside the plotting region 
# inset = c(x,y): how to move the legend relative to the 'bottom' location

#####################################################################################
######################### FIGURE 5: FINAL POPULATION SIZES ##########################
#####################################################################################

## Select file name:
file_name <- "Fig5_FinalPopSize_LM_IBM.pdf"
setwd(workingDir_Figures)

pdf(file_name, width=12, height=6)
layout(matrix(c(1,2,3,4,5,6,7,7,7), ncol=3, byrow=TRUE), heights=c(0.45,0.45,0.1))
par(oma=c(1, 1, 1, 1), mar=c(4, 4, 2, 2))

t <- boxplot(finalPop3$A, finalPop3$B, finalPop3$C, finalPop3$D, finalPop3$E, finalPop3$F, finalPop3$G, finalPop3$H, finalPop3$I, names=LETTERS[1:9], pch=16, cex=1.2, col=colMTN1, outcol=colMTN1_alpha, varwidth=FALSE, medlwd=1, medcol="white", boxlty=1, whisklty=1, staplelty=0, ylim=c(0,200), yaxt = "n", axes=FALSE)
box(bty="l")
axis(2, font.lab=2, at=seq(0, 200, by=50), labels=seq(0, 200, by=50), las=2)
axis(1, font.lab=2, at=seq(1, 9, by=1), labels=LETTERS[1:9], las=1)
title(main="Mountain gorilla, 1%", font.lab=2)

t <- boxplot(finalPop2$A, finalPop2$B, finalPop2$C, finalPop2$D, finalPop2$E, finalPop2$F, finalPop2$G, finalPop2$H, finalPop2$I, names=LETTERS[1:9], pch=18, cex=1.2, col=colMTN2, outcol=colMTN2_alpha, varwidth=FALSE, medlwd=1, medcol="white", boxlty=1, whisklty=1, staplelty=0, ylim=c(0,200), yaxt = "n", axes=FALSE)
box(bty="l")
axis(2, font.lab=2, at=seq(0, 200, by=50), labels=seq(0, 200, by=50), las=2)
axis(1, font.lab=2, at=seq(1, 9, by=1), labels=LETTERS[1:9], las=1)
title(main="Mountain gorilla, 2%", font.lab=2)

t <- boxplot(finalPop1$A, finalPop1$B, finalPop1$C, finalPop1$D, finalPop1$E, finalPop1$F, finalPop1$G, finalPop1$H, finalPop1$I, names=LETTERS[1:9], pch=17, cex=1.2, col=colMTN3, outcol=colMTN3_alpha, varwidth=FALSE, medlwd=1, medcol="white", boxlty=1, whisklty=1, staplelty=0, ylim=c(0,200), yaxt = "n", axes=FALSE)
box(bty="l")
axis(2, font.lab=2, at=seq(0, 200, by=50), labels=seq(0, 200, by=50), las=2)
axis(1, font.lab=2, at=seq(1, 9, by=1), labels=LETTERS[1:9], las=1)
title(main="Mountain gorilla, 3.2%", font.lab=2)

t <- boxplot(finalPop7$A, finalPop7$B, finalPop7$C, finalPop7$D, finalPop7$E, finalPop7$F, finalPop7$G, finalPop7$H, finalPop7$I, names=LETTERS[1:9], pch=16, cex=1.2, col=colMTN1, outcol=colMTN1_alpha, varwidth=FALSE, medlwd=1, medcol="white", boxlty=1, whisklty=2, staplelty=0, ylim=c(0,200), yaxt = "n", axes=FALSE)
box(bty="l")
axis(2, font.lab=2, at=seq(0, 200, by=50), labels=seq(0, 200, by=50), las=2)
axis(1, font.lab=2, at=seq(1, 9, by=1), labels=LETTERS[1:9], las=1)

t <- boxplot(finalPop6$A, finalPop6$B, finalPop6$C, finalPop6$D, finalPop6$E, finalPop6$F, finalPop6$G, finalPop6$H, finalPop6$I, names=LETTERS[1:9], pch=18, cex=1.2, col=colMTN2, outcol=colMTN2_alpha, varwidth=FALSE, medlwd=1, medcol="white", boxlty=1, whisklty=2, staplelty=0, ylim=c(0,200), yaxt = "n", axes=FALSE)
box(bty="l")
axis(2, font.lab=2, at=seq(0, 200, by=50), labels=seq(0, 200, by=50), las=2)
axis(1, font.lab=2, at=seq(1, 9, by=1), labels=LETTERS[1:9], las=1)

t <- boxplot(finalPop5$A, finalPop5$B, finalPop5$C, finalPop5$D, finalPop5$E, finalPop5$F, finalPop5$G, finalPop5$H, finalPop5$I, names=LETTERS[1:9], pch=17, cex=1.2, col=colMTN3, outcol=colMTN3_alpha, varwidth=FALSE, medlwd=1, medcol="white", boxlty=1, whisklty=2, staplelty=0, ylim=c(0,200), yaxt = "n", axes=FALSE)
box(bty="l")
axis(2, font.lab=2, at=seq(0, 200, by=50), labels=seq(0, 200, by=50), las=2)
axis(1, font.lab=2, at=seq(1, 9, by=1), labels=LETTERS[1:9], las=1)

mtext("Reintroduction Scenario", side=1, line=-5, outer=TRUE, cex=1, font=2)
mtext("Population size after 50 years", side=2, line=-1, outer=TRUE, cex=1, font=2, las=0)

par(mai=c(0,0,0,0))
plot.new()
legend("bottom", legend=c("Leslie Matrix", "Individual-based Model"), lty=c(1,2), text.font=2, xpd = TRUE, horiz = FALSE, inset = c(0, 0.02), bty = "y")
dev.off()

file_name <- "Fig5_FinalPopSize_LM_IBM.pdf"
setwd(workingDir_Figures)

pdf(file_name, width=12, height=6)
layout(matrix(c(1,2,3,4,5,6,7,7,7), ncol=3, byrow=TRUE), heights=c(0.4,0.4,0.2))
par(oma=c(1, 1, 1, 1), mar=c(4, 4, 2, 2))

boxplot(finalPop7$A, finalPop3$A, 
        finalPop7$B, finalPop3$B, 
        finalPop7$C, finalPop3$C, 
        finalPop7$D, finalPop3$D, 
        finalPop7$E, finalPop3$E, 
        finalPop7$F, finalPop3$F, 
        finalPop7$G, finalPop3$G, 
        finalPop7$H, finalPop3$H, 
        finalPop7$I, finalPop3$I, 
        pch=15, col=colMTN1, outcol=colMTN1_alpha, 
        medlwd=1, medcol="white", boxlty=1, whisklty=1, staplelty=0, 
        varwidth=FALSE, ylim=c(0,100), yaxt = "n")
rect(1.5,-1,2.5,220, col=colWLG_alpha,lty=0)
rect(3.5,-1,4.5,220, col=colWLG_alpha,lty=0)
rect(5.5,-1,6.5,220, col=colWLG_alpha,lty=0)
rect(7.5,-1,8.5,220, col=colWLG_alpha,lty=0)
rect(9.5,-1,10.5,220, col=colWLG_alpha,lty=0)
rect(11.5,-1,12.5,220, col=colWLG_alpha,lty=0)
rect(13.5,-1,14.5,220, col=colWLG_alpha,lty=0)
rect(15.5,-1,16.5,220, col=colWLG_alpha,lty=0)
rect(17.5,-1,18.5,220, col=colWLG_alpha,lty=0)
boxplot(finalPop7$A, finalPop3$A, 
        finalPop7$B, finalPop3$B, 
        finalPop7$C, finalPop3$C, 
        finalPop7$D, finalPop3$D, 
        finalPop7$E, finalPop3$E, 
        finalPop7$F, finalPop3$F, 
        finalPop7$G, finalPop3$G, 
        finalPop7$H, finalPop3$H, 
        finalPop7$I, finalPop3$I, 
        pch=15, col=colMTN1, outcol=colMTN1_alpha, 
        medlwd=1, medcol="white", boxlty=1, whisklty=1, staplelty=0, 
        varwidth=FALSE, ylim=c(0,100), yaxt = "n", add=TRUE)

boxplot(finalPop3$A, finalPop3$B, finalPop3$C, finalPop3$D, finalPop3$E, finalPop3$F, finalPop3$G, finalPop3$H, finalPop3$I, names=LETTERS[1:9], pch=15, col=colMTN1_alpha, outcol=colMTN1_alpha, varwidth=FALSE, medlwd=1, medcol="white", boxlty=1, whisklty=1, staplelty=0, ylim=c(0,220), yaxt = "n", add=TRUE)
axis(2, font.lab=2, at=seq(0, 220, by=50), labels=seq(0, 200, by=50), las=2)
lines(x=c(-5:51), y=rep(50, 57), col=colN50, lwd=1, lty=1)
title(main="Mountain gorilla, 1%", xlab="Reintroduction Scenario", ylab="Population size after 50 years", font.lab=2)

boxplot(finalPop2$A, finalPop2$B, finalPop2$C, finalPop2$D, finalPop2$E, finalPop2$F, finalPop2$G, finalPop2$H, finalPop2$I, names=LETTERS[1:9], pch=18, col=colMTN2, outcol=colMTN2_alpha, varwidth=FALSE, medlwd=1, medcol="white", boxlty=1, whisklty=1, staplelty=0, ylim=c(0,220), yaxt = "n")
axis(2, font.lab=2, at=seq(0, 220, by=50), labels=seq(0, 200, by=50), las=2)
lines(x=c(-5:51), y=rep(50, 57), col=colN50, lwd=1, lty=1)
title(main="Mountain gorilla, 2%", xlab="Reintroduction Scenario", ylab="Population size after 50 years", font.lab=2)

boxplot(finalPop1$A, finalPop1$B, finalPop1$C, finalPop1$D, finalPop1$E, finalPop1$F, finalPop1$G, finalPop1$H, finalPop1$I, names=LETTERS[1:9], pch=17, col=colMTN3, outcol=colMTN3_alpha, varwidth=FALSE, medlwd=1, medcol="white", boxlty=1, whisklty=1, staplelty=0, ylim=c(0,220), yaxt = "n")
axis(2, font.lab=2, at=seq(0, 220, by=50), labels=seq(0, 200, by=50), las=2)
lines(x=c(-5:51), y=rep(50, 57), col=colN50, lwd=1, lty=1)
title(main="Mountain gorilla, 3.2%", xlab="Reintroduction Scenario", ylab="Population size after 50 years", font.lab=2)

boxplot(finalPop6$A, finalPop6$B, finalPop6$C, finalPop6$D, finalPop6$E, finalPop6$F, finalPop6$G, finalPop6$H, finalPop6$I, names=LETTERS[1:9], pch=18, col=colMTN2, outcol=colMTN2_alpha, varwidth=FALSE, medlwd=1, medcol="white", boxlty=1, whisklty=2, staplelty=0, ylim=c(0,220), yaxt = "n")
axis(2, font.lab=2, at=seq(0, 220, by=50), labels=seq(0, 200, by=50), las=2)
lines(x=c(-5:51), y=rep(50, 57), col=colN50, lwd=1, lty=1)
title(xlab="Reintroduction Scenario", ylab="Population size after 50 years", font.lab=2)

boxplot(finalPop5$A, finalPop5$B, finalPop5$C, finalPop5$D, finalPop5$E, finalPop5$F, finalPop5$G, finalPop5$H, finalPop5$I, names=LETTERS[1:9], pch=17, col=colMTN3, outcol=colMTN3_alpha, varwidth=FALSE, medlwd=1, medcol="white", boxlty=1, whisklty=2, staplelty=0, ylim=c(0,220), yaxt = "n")
axis(2, font.lab=2, at=seq(0, 220, by=50), labels=seq(0, 200, by=50), las=2)
lines(x=c(-5:51), y=rep(50, 57), col=colN50, lwd=1, lty=1)
title(xlab="Reintroduction Scenario", ylab="Population size after 50 years", font.lab=2)

par(mai=c(0,0,0,0))
plot.new()
legend("center", legend=c("Mountain Gorillas - 1%", "Mountain Gorillas - 2%", "Mountain Gorillas - 3.2%"), fill=c(colMTN1, colMTN2, colMTN3), text.font=2, cex=1, xpd = TRUE, horiz = FALSE)
dev.off()

#####################################################################################
############ SUPPLEMENTARY FIGURES 1-4: LESLIE MATRIX PROJECTION PLOTS ##############
#####################################################################################

## Select file name:
file_name <- "SuppFig3_LM_Projections_MTN3.pdf"
file_name <- "SuppFig2_LM_Projections_MTN2.pdf"
file_name <- "SuppFig1_LM_Projections_MTN1.pdf"
file_name <- "SuppFig_LM_Projections_WLG.pdf"

setwd(workingDir_Output)
## Select exticntion risk files:
probExt_lm <- read.csv("Results/Results_LM_mtn_3%.csv")
probExt_lm <- read.csv("Results/Results_LM_mtn_2%.csv")
probExt_lm <- read.csv("Results/Results_LM_mtn_1%.csv")
probExt_lm <- read.csv("Results/Results_LM_wlg.csv")

## Select the correct folder for either WLG or MTN data
workingDir_LM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/PVA_Output/LM_Projection_50year_mtn_3%/")
workingDir_LM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/PVA_Output/LM_Projection_50year_mtn_2%/")
workingDir_LM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/PVA_Output/LM_Projection_50year_mtn_1%/")
workingDir_LM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/PVA_Output/LM_Projection_50year_wlg/")

setwd(workingDir_LM)
allScenarioFiles <- list.files(pattern="*.csv")

for (i in 1:length(allScenarioFiles)){
  assign(allScenarioFiles[i], 
         read.csv(allScenarioFiles[i], header=TRUE)
  )}

stochObjects <- c("LM_Stoch_Scenario1.csv","LM_Stoch_Scenario2.csv","LM_Stoch_Scenario3.csv","LM_Stoch_Scenario4.csv","LM_Stoch_Scenario5.csv","LM_Stoch_Scenario6.csv","LM_Stoch_Scenario7.csv","LM_Stoch_Scenario8.csv","LM_Stoch_Scenario9.csv")
detObjects <- c("LM_Det_Scenario1.csv","LM_Det_Scenario2.csv","LM_Det_Scenario3.csv","LM_Det_Scenario4.csv", "LM_Det_Scenario5.csv","LM_Det_Scenario6.csv","LM_Det_Scenario7.csv","LM_Det_Scenario8.csv","LM_Det_Scenario9.csv")

setwd(workingDir_Figures)
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
  lines(apply(tempX, 1, mean)~time, type="l", col="white", lwd=4) ## plot mean projection from stochastic LM simulations
  lines(N_projected_detX~time, type="l", col=1, lwd=2, lty=4) ## replot deterministic projection
  title(main=paste0("Scenario ",scenario), sub=paste0("Extinction Risk = ",probExt, "%"), cex.main=1, cex.sub=1, col.sub=1, font.sub=3)
  qtiles <- apply(tempX, 1, function(v) quantile(v, probs=c(0.05, 0.95))) ## plot 95% confidence intervals for simulations
  lines((0:(nrow(tempX)-1)), qtiles[1,], col=colMTN1, lwd=3, lty=1)
  lines((0:(nrow(tempX)-1)), qtiles[2,], col=colMTN1, lwd=3, lty=1)
  lines(x=c(-5:50), y=rep(50, 56), col=colN50, lwd=1, lty=1) ## add a line for the 50 individual mark
}
dev.off()

#####################################################################################
################## SUPPLEMENTARY FIGURE 5-8: IBM PROJECTION PLOTS ###################
#####################################################################################

setwd(workingDir_Output)
## Select file name:
file_name <- "SuppFig6_IBM_Projections_MTN3.pdf"
file_name <- "SuppFig5_IBM_Projections_MTN2.pdf"
file_name <- "SuppFig4_IBM_Projections_MTN1.pdf"
file_name <- "SuppFig_IBM_Projections_WLG.pdf"

## Select exticntion risk files:
probExt_ibm <- read.csv("Results/Results_IBM_mtn_3%.csv")
probExt_ibm <- read.csv("Results/Results_IBM_mtn_2%.csv")
probExt_ibm <- read.csv("Results/Results_IBM_mtn_1%.csv")
probExt_ibm <- read.csv("Results/Results_IBM_wlg.csv")

## Select the correct folder for either WLG or MTN data
workingDir_IBM <- ("IBM_Projection_50year_mtn_3%")
workingDir_IBM <- ("IBM_Projection_50year_mtn_2%")
workingDir_IBM <- ("IBM_Projection_50year_mtn_1%")
workingDir_IBM <- ("IBM_Projection_50year_wlg")

setwd(workingDir_IBM)
allScenarioFiles <- list.files(pattern="*.csv")

for (i in 1:length(allScenarioFiles)){
  assign(allScenarioFiles[i], 
         read.csv(allScenarioFiles[i], header=TRUE)
  )}

stochObjects <- c("IBM_Scenario1.csv","IBM_Scenario2.csv","IBM_Scenario3.csv","IBM_Scenario4.csv","IBM_Scenario5.csv","IBM_Scenario6.csv","IBM_Scenario7.csv","IBM_Scenario8.csv","IBM_Scenario9.csv")

timeunit <- 1/12 ## time interval for the plots

setwd(workingDir_Figures)
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
  lines((0:(nrow(resX)-1))*timeunit, apply(resX, 1, mean), type="l", col="white", lwd=3)
  
  ## add 95% upper/lower limits
  qtiles <- apply(resX, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
  lines((0:(nrow(resX)-1))*timeunit, qtiles[1,], col=colMTN1, lwd=1, lty=1)
  lines((0:(nrow(resX)-1))*timeunit, qtiles[2,], col=colMTN1, lwd=1, lty=1)
  
  title(main=paste0("Scenario ",scenario), sub=paste0("Extinction Risk = ",probExt, "%"), cex.main=1, cex.sub=1, col.sub=1,  font.sub=3)
  
  ## add a line for the 50 individual mark
  lines(x=c(-5:50), y=rep(50, 56), col=colN50, lwd=1, lty=1)
}
dev.off()
