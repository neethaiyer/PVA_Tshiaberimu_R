#######################################################
############### DETERMINE CORRECT alpha ###############
#######################################################

## Make sure your current working directory is the main folder:
setwd("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/")

## If these directories do not exsit, run this:
dir.create("PVA_Output")
dir.create("PVA_Input")
dir.create("PVA_Output/Results")
dir.create("PVA_Output/LM_Projection")
dir.create("PVA_Output/IBM_Projection")

## Select the working directory for input files:
workingDir_Input <- "/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/PVA_Input/"

## Select the working directory for output files:
workingDir_Output <- "/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/PVA_Output/"

## DAMIEN !!!!! Your wd: ## Make sure to update as above!
## workingDir <- "~/Documents/git repositories/PVA_Tshiaberimu_R/"

## Source the functions used in simulations below:
source("1. Function Definitions.R") 

####################################################################################
############################ READ LIFE HISTORY TABLE ###############################
####################################################################################
setwd(workingDir_Input)

## Your life history tables should have at least 3 columns: age, mortality rate, and fertility rate:
dat <- read.csv("Gorilla_LifeTables.csv")
dat$fertilityrate_2percent <- dat[,3]*.786 
## fertility rates multiplied by factor less than 1 to get eigen values of 1.01 which corresponds to a 1% growth rate
dat$fertilityrate_1percent <- dat[,3]*.643 
## fertility rates multiplied by factor less than 1 to get eigen values of 1.02 which corresponds to a 2% growth rate

##################################################################
############## Create an object that selects the LM ##############
##################################################################

selectLM <- read.csv("LeslieMatrix_WLG.csv")

###############################################################################
############## SET THE INITIAL CONDITIONS OF THE LM & IBM MODELS ##############
###############################################################################

## Reintroduction Scenarios:
ReintroScenario <- read.csv("ReintroductionScenarios_LM.csv") ## csv file with Reintroduction Scenarios for LM

## Leslie Matrix parameters:
mat <- as.matrix(selectLM) ## LM needs to be converted to matrix object!

## Time parameters:
nyears <- 100 ## Projection Period
nruns <- 1000 ## Number of simulations to run
timeunit <- 1/12 ## timestep for IBM

## Initial demographic parameters: survivorship, fertility, and weaning age
datX <- dat[,c(1,4:5)] ## Subset appropriate life history columns: dat[,c(1,4:5)] for WLG, dat[,1:3] for MTN
## NOTE: this subsetting is needed because columns for dat are specified in FUNCTIONS 8 and 9
weaningAge <- 4.5 ## 4.5 for WLG, 3.5 for MTN
adultAge <- 10 ## 10 for WLG, 8 for MTN

temp <- pop_projection(tfinal=nyears, LM=mat, No=ReintroScenario[,11]) ## run the LM projection
Nfinal <- data.frame(ReintroScenario$age, temp[,51])
colnames(Nfinal) <- c("age","N")
sum(Nfinal$N) ## check if Nfinal is about the same as Nfinal_expected for the required lambda value
Nfinal_expected <- 1.032^50*30 ## Nfinal_expected = lambda^t*No (lambda for 3% = 1.032)
N_random <- rep(Nfinal$age, round(Nfinal$N)) ## create vector with individuals and their ages
N_random <- sample(N_random, size=100, replace = FALSE, prob = NULL) ## randomly sample 100 individuals in each age class
N_random <- data.frame(subset(N_random, N_random>weaningAge)) ## remove individuals younger than weaning age
colnames(N_random) <- "TestScenario"

growthRates <- data.frame(alpha_value = NA, 
                          growth_rate = NA, lambda = NA)

## run the IBM projection
timeunit<-1/12
initalConditions <- convertToList(scenario = N_random, adultAge=adultAge, weaningAge=weaningAge) ## define initial conditions based on ages of females randomly sampled earlier in N_random
nruns <- 10
alpha <- 0.10 ## set alpha value

res <- matrix(0, nrow=trunc(nyears/timeunit)+1, ncol=nruns)
for(j in 1:length(initalConditions)){
  for(i in 1:nruns){
    print(i)
    abmDataLog <- simTshia(ages0 = initalConditions[[j]][,1], status0 = initalConditions[[j]][,2], time0 = initalConditions[[j]][,3], nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F)
    nindiv <- tapply(abmDataLog$status,abmDataLog$timestep, function(v) length(v)+rbinom(1, sum(v=="L"), .5))##we're adding the unweaned females
    res[1:length(nindiv),i] <- nindiv
  }
}

finalPop <- as.numeric(res[nrow(res),])
startPop <- as.numeric(res[840,])
##startPop <- rep(100,nruns) ## if using No at t=0

logLambda <- mean((1/30)*log(finalPop/startPop)) ## nyears for the census time period, loglambda = 1/timeperiod*log(Ntfinal)/Nt0
lambda <- exp(logLambda)
growthRates[2,1:3] <- c(alpha, round(logLambda, digits=3), round(lambda, digits=3))
growthRates <- growthRates[order(-growthRates$alpha_value),] 
growthRates

setwd("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/")
write.csv(growthRates, file="growthRates_LM2%_100yrs_FINAL.csv", row.names=F)

## Examine how alpha and growth rate vary:
growthRates1 <- read.csv("growthRates_LM2%_100yrs_FINAL.csv")
growthRates2 <- read.csv("growthRates_LM2%.csv")
growthRates3 <- read.csv("growthRates_LM3%.csv")

plot(growthRates3[,1:2], col=2)
points(growthRates2[,1:2], col=4)
points(growthRates1[,1:2])
abline(lm(growthRates3[,2]~growthRates3[,1]), col=2)
abline(lm(growthRates2[,2]~growthRates2[,1]), col=4)
abline(lm(growthRates1[,2]~growthRates1[,1]))




