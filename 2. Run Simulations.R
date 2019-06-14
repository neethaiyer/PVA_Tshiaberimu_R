## Set the working directory:
workingDir <- "~/Box Sync/PVA_Paper/PVA_Tshiaberimu_R/"
##workingDir <- "~/Documents/git repositories/PVA_Tshiaberimu_R/"

source("1. Function Definition.R")

#####################################################################################
############ PART 1: Historical population trends for Tshiaberimu gorillas ##########
#####################################################################################

## Take a look at the historical population trajectories using Tshiaberimu census data
year <- c(1959,1986,1995,1996,2003,2004,2006,2007,2008,2009,2011,2012,2013,2016,2017) ## census years
numyears <- 1959:2017 ## census time period
N <- c(35,20,17,16,20,20,21,22,18,16,6,6,7,6,6) ## census data

## let's look at the rate of change in this population
## lambda is the finite rate of increase of a population over one time step. r is the intrinsinc rate of growth. negative r values indicate a population in decline. lambda < 1 indicates a decline. the relationship between lambda and r : lambda = Nt+1  / Nt, r = ln(lambda), lambda = e^r
logLambda <- (1/58)*log(6/35) ## 58 years for the census time period, lambda = 1/timeperiod*log(Ntfinal)/Nt0
lambda <- exp(logLambda)
popEst <- 35*(exp(logLambda))^(0:58) ## this is the expected rate of change in the population given Ntfinal and Nt0
popEst ## these are the predicted population estimates given the calculated lambda value

## let's fit these parameter estimates to a linear model to calculate the r and lambda values to get a more accurate estimate of these parameters:
modelGeom <- lm(log(N)~year) ## should be linear on a log scale
r_lm <- modelGeom$coef[2] ## take the slope of the line from this linear model for the intrinsic rate of growth r=-0.0289175
lambda_lm <- exp(modelGeom$coef[2]) ## lambda=0.9714966
popEst_lm <- 35*(exp(r_lm))^(0:58)

## plot the actual population sizes from census data and the expected population size:
plot(year, N, xlab="Census Year", 
     ylab="Estimated population size, N", 
     pch=19, type="o",
     ylim=c(0,50), 
     xlim=c(numyears[1],numyears[length(numyears)]), 
     las=1, cex.main=0.8, cex.lab=0.8, cex.axis=0.8, font.lab=2)
lines(numyears, popEst_lm, col=2, lty=2, lwd=2)

#####################################################################################
##################### PART 2: Leslie Matrix model (Simple PVA) ######################
#####################################################################################

## So, the population is clearly declining. Let's estimate population trajectories based on reintroduction scenarios via Leslie Matrix estimation:

## Data from Breuer et al (2010), Breuer (2008) for western lowland gorillas (WLG)
## Note the first age of reproduction is 10 years and the fertility rate was constant at 0.16 (16%) for all reproductive females. 

dat <- read.csv(paste0(workingDir, "LifeTable_WLG_Breuer.csv"))

## Data from Bronikowski et al (2016) for mountain gorillas (MTN)
## Note the first age of reproduction is 8 years although the first birth is usually at 10 years old for mountain gorillas and the fertility rate varied for each adult year.
dat1 <- read.csv(paste0(workingDir, "LifeTable_MTN_Bronikowski.csv"))
dat1$fertilityrate_1percent <- dat1[,3]*.643 ## fertility rates multiplied by factor less than 1 to get eigen values of 1.01 which corresponds to a 1% growth rate
dat1$fertilityrate_2percent <- dat1[,3]*.786 ## fertility rates multiplied by factor less than 1 to get eigen values of 1.02 which corresponds to a 2% growth rate