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

## Projections were conducted using 4 possible growth rates. Make sure the right csv file is selected for each simulation. 
## A. Mountain gorillas with fertility rates that correspond to 3% growth rate
## B. Mountain gorillas with fertility rates that correspond to 2% growth rate
## C. Mountain gorillas with fertility rates that correspond to 1% growth rate
## D. Western lowland gorillas with constant fertility rates (due to limited data)
## Age of first reproduction = 8 years (for MTN) and 10 years (for WLG)

## Mountain gorilla (MTN) life history parameters taken from Bronikowski et al (2016)
## Western lowland gorilla (WLG) life history parameters taken from Breuer et al (2010), Breuer (2008) 

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

####################################################################################
################# OPTIONAL: CREATE CSV FILES with LESLIE MATRICES ##################
####################################################################################

## If you don't already have Leslie Matrices, create using the function leslieMatrix. Make sure you select the correct columns in your Life Table file (you need: age, mortality rate, and fertility rate):
leslieMatrix(lifetable=dat[,1:3], filename="LeslieMatrix_MTN_3%.csv")
leslieMatrix(lifetable=dat[,c(1:2, 6)], filename="LeslieMatrix_MTN_2%.csv")
leslieMatrix(lifetable=dat[,c(1:2, 7)], filename="LeslieMatrix_MTN_1%.csv")
leslieMatrix(lifetable=dat[,c(1, 4:5)], filename="LeslieMatrix_WLG.csv")

##################################################################
############## Create an object that selects the LM ##############
##################################################################

selectLM <- read.csv("LeslieMatrix_MTN_1%.csv")

###############################################################################
############## SET THE INITIAL CONDITIONS OF THE LM & IBM MODELS ##############
###############################################################################

## Reintroduction Scenarios:
ReintroScenario <- read.csv("ReintroductionScenarios_LM.csv") ## csv file with Reintroduction Scenarios for LM
ReintroScenario_IBM <- read.csv("ReintroductionScenarios_IBM.csv") ## csv file with Reintroduction Scenarios for IBM

## Leslie Matrix parameters:
mat <- as.matrix(selectLM) ## LM needs to be converted to matrix object!

## Time parameters:
nyears <- 50 ## Projection Period
nruns <- 1000 ## Number of simulations to run
timeunit <- 1/12 ## timestep for IBM

## Initial demographic parameters: survivorship, fertility, and weaning age
datX <- dat[,1:3] ## Subset appropriate life history columns: dat[,c(1,4:5)] for WLG, dat[,1:3] for MTN
## NOTE: this subsetting is needed because columns for dat are specified in FUNCTIONS 8 and 9
weaningAge <- 3.5 ## 4.5 for WLG, 3.5 for MTN
adultAge <- 8 ## 10 for WLG, 8 for MTN
alpha <- 0.60 ## see "calculate_alpha_value.R" for more details

## Depending on the adult female age and weaning age, create a list with the starting conditions for each scenario of the IBM
initalConditions <- convertToList(scenario = ReintroScenario_IBM, adultAge=adultAge, weaningAge=weaningAge)

#####################################################################################
##################### PART 1: Leslie Matrix model (Simple PVA) ######################
#####################################################################################
setwd(workingDir_Output)

###############################################################################
######################### RUN THE LESLIE MATRIX MODELS ########################
###############################################################################

## Apply LM projection functions using each reintrodcution scenario
## First, use the deterministic function:
  projectPop <- for(i in 2:ncol(ReintroScenario)){
    No <- ReintroScenario[,i] ## Get the reintroduction scenario
    N <- pop_projection(tfinal=nyears, LM=mat, No=No) ## Apply pop_projection function to No for this scenario
    scenario <- strsplit(colnames(ReintroScenario)[i], "_")[[1]][2] ## Get the last element of the column name for each reintroducion scenario
    det <- assign(paste0("N_projected_det", scenario), apply(N,2,sum))  ## The assign function takes a variable name as a character string and assigns a value to it. In this case, the values are N at each time step of the projection
    write.csv(det, file=paste0("LM_Projection/LM_Det_Scenario", i-1, ".csv"), row.names=F)
  }

## Second, use the stochastic function:
## create an empty matrix that will save the number of individuals for each year of the projection for each run of the LM projection 
temp <- matrix(0, nrow=nyears+1, ncol=nruns)

for(j in 1:(length(ReintroScenario)-1)){
  for(i in 1:nruns) {
    temp[1:(nyears+1),i] <- apply(stoch_projection(tfinal=nyears, LM=mat, No=ReintroScenario[,j+1]),2,sum)
  }
  write.csv(temp, file=paste0("LM_Projection/LM_Stoch_Scenario", j,".csv"), row.names=F)
}

###############################################################################
############## CALCULATE RISK OF EXTINCTION & Likelihood Pop = 50 #############
###############################################################################
### Now that the csv files have been written, we may not want to re-run the code for as long in the future, so we can just read the generated files and plot the data directly. Make sure you're reading the csv files from the correct folder. You might want to move the files from the LM_Projection folder into a new folder. 

setwd(workingDir_Output)
## Select the correct folder:
##workingDir_LM <- "LM_Projection_50year_mtn_3%"
##workingDir_LM <- "LM_Projection_50year_mtn_2%"
workingDir_LM <- "LM_Projection_50year_mtn_1%"
##workingDir_LM <- "LM_Projection_50year_wlg"

setwd(workingDir_LM)
allScenarioFiles <- list.files(pattern="LM_Stoch.*\\.csv")

## Calculate the extinction risk at the end of 50 years for all nruns of simulation.
## Calculate the likelihood that the population reaches at least 50 individuals within 50 years.

results_LM <- data.frame(scenario = as.factor(LETTERS[1:length(allScenarioFiles)]), 
                         prob_50 = NA, extn_Risk = NA)
index <- 0
for(i in 1:length(allScenarioFiles)){
  index <- index+1
  tempx <- as.matrix(read.csv(allScenarioFiles[[i]]))
  ext <- tempx[nrow(tempx),]==0
  extRisk <- mean(ext)
  probNe_50 <- mean(tempx[nrow(tempx),]>=50)
  results_LM[index,2:3] <- c(probNe_50, extRisk)*100	
}

setwd(workingDir_Output)
## Write csv files to save the results:
##write.csv(results_LM, file="Results/Results_LM_MTN_3%.csv", row.names=F)
##write.csv(results_LM, file="Results/Results_LM_MTN_2%.csv", row.names=F)
write.csv(results_LM, file="Results/Results_LM_MTN_1%.csv", row.names=F)
##write.csv(results_LM, file="Results/Results_LM_WLG.csv", row.names=F)

#####################################################################################
################ PART 2: Individual Based Model (IBM) (Complex PVA) #################
#####################################################################################
setwd(workingDir_Output)
###############################################################################
################################## RUN THE IBM ################################
###############################################################################

res <- matrix(0, nrow=trunc(nyears/timeunit)+1, ncol=nruns)
for(j in 1:length(initalConditions)){
  for(i in 1:nruns){
    print(i)
    abmDataLog <- simTshia(ages0 = initalConditions[[j]][,1], status0 = initalConditions[[j]][,2], time0 = initalConditions[[j]][,3], nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F)
    nindiv <- tapply(abmDataLog$status,abmDataLog$timestep, function(v) length(v)+rbinom(1, sum(v=="L"), .5))##we're adding the unweaned females
    res[1:length(nindiv),i] <- nindiv
  }
  write.csv(res, file=paste0("IBM_Projection/IBM_Scenario", j,".csv"), row.names=F)
}

###############################################################################
################ CALCULATE RISK OF EXTINCTION & Likelihood Pop = 50 ###########
###############################################################################
### Now that the csv files have been written, we may not want to re-run the code for as long in the future, so we can just read the generated files and plot the data directly. Make sure you're reading the csv files from the correct folder. You might want to move the files from the IBM_Projection folder into a new folder. 

setwd(workingDir_Output)
## Select the correct folder:
workingDir_IBM <- "IBM_Projection_50year_mtn_3%"
##workingDir_IBM <- "IBM_Projection_50year_mtn_2%"
##workingDir_IBM <- "IBM_Projection_50year_mtn_1%"
##workingDir_IBM <- "IBM_Projection_50year_wlg"

setwd(workingDir_IBM)
allScenarioFiles <- list.files(pattern="*.csv")

## Calculate the extinction risk at the end of 50 years for all nruns of simulation.
## Calculate the likelihood that the population reaches at least 50 individuals within 50 years.

results_IBM <- data.frame(scenario = as.factor(LETTERS[1:length(allScenarioFiles)]), 
                           probNe_50 = NA, 
                           extn_Risk = NA)

index <- 0
for(i in 1:length(allScenarioFiles)){
  index <- index+1
  res <- as.matrix(read.csv(paste0(workingDir_IBM, allScenarioFiles[[i]])))
  probNe_50 <- mean(res[nrow(res),]>=50)
  extn_Risk <- mean(res[nrow(res),]==0)
  results_IBM[index,2:3] <- c(probNe_50, extn_Risk)*100
}

setwd(workingDir_Output)
## Write csv files to save the results:
write.csv(results_IBM, file="Results/extn_ibm_MTN_3%.csv", row.names=F)
#write.csv(results_IBM, file="Results/extn_ibm_MTN_2%.csv", row.names=F)
#write.csv(results_IBM, file="Results/extn_ibm_MTN_1%.csv", row.names=F)
#write.csv(results_IBM, file="Results/extn_ibm_WLG.csv", row.names=F)

###############################################################################
################# CREATE CSV FILES FOR FINAL POPULATION SIZES #################
###############################################################################
setwd(workingDir_Output)
## Select the correct folder for either LM or IBM:
workingDir_LM <- "LM_Projection_50year_mtn_3%"
##workingDir_LM <- "LM_Projection_50year_mtn_2%"
##workingDir_LM <- "LM_Projection_50year_mtn_1%"
##workingDir_LM <- "LM_Projection_50year_wlg"
workingDir_IBM <- "IBM_Projection_50year_mtn_3%"
##workingDir_IBM <- "IBM_Projection_50year_mtn_2%"
##workingDir_IBM <- "IBM_Projection_50year_mtn_1%"
##workingDir_IBM <- "IBM_Projection_50year_wlg"

setwd(workingDir_LM)
##setwd(workingDir_IBM)

allScenarioFiles <- list.files(pattern="*.csv")
for (i in 1:length(allScenarioFiles)){
  assign(allScenarioFiles[i], 
         read.csv(allScenarioFiles[i], header=TRUE)
  )}

########## IMPORTANT !!!!!! #########
########## IMPORTANT !!!!!! #########
########## IMPORTANT !!!!!! #########
## Select simulation objects from list. Uncommented depending on whether IBM or LM:
stochObjects <- c("LM_Stoch_Scenario1.csv","LM_Stoch_Scenario2.csv","LM_Stoch_Scenario3.csv","LM_Stoch_Scenario4.csv","LM_Stoch_Scenario5.csv","LM_Stoch_Scenario6.csv","LM_Stoch_Scenario7.csv","LM_Stoch_Scenario8.csv","LM_Stoch_Scenario9.csv")
##stochObjects <- c("IBM_Scenario1.csv","IBM_Scenario2.csv","IBM_Scenario3.csv","IBM_Scenario4.csv","IBM_Scenario5.csv","IBM_Scenario6.csv","IBM_Scenario7.csv","IBM_Scenario8.csv","IBM_Scenario9.csv")

Nfinal <- data.frame(matrix(ncol=length(stochObjects), nrow=nruns))
colnames(Nfinal) <- LETTERS[1:length(stochObjects)]

index <- 0
for(j in 1:length(stochObjects)){
  index<-index+1
  resX <- get(stochObjects[j])
  Nfinal[,j] <- as.numeric(resX[nrow(resX),])
}

##write.csv(Nfinal, file="Results_LM_Nfinal_MTN_3%.csv", row.names=F)
##write.csv(Nfinal, file="Results_LM_Nfinal_MTN_2%.csv", row.names=F)
##write.csv(Nfinal, file="Results_LM_Nfinal_MTN_1%.csv", row.names=F)
##write.csv(Nfinal, file="Results_LM_Nfinal_WLG.csv", row.names=F)
##write.csv(Nfinal, file="Results_IBM_Nfinal_MTN_3%.csv", row.names=F)
##write.csv(Nfinal, file="Results_IBM_Nfinal_MTN_2%.csv", row.names=F)
##write.csv(Nfinal, file="Results_IBM_Nfinal_MTN_1%.csv", row.names=F)
##write.csv(Nfinal, file="Results_IBM_Nfinal_WLG.csv", row.names=F)

###############################################################################
################### CALCULATE GROWTH RATE R FOR EACH MODEL ####################
###############################################################################
setwd(workingDir_Output)
## Select the correct folder for either LM or IBM:
workingDir_LM <- "LM_Projection_50year_mtn_3%"
##workingDir_LM <- "LM_Projection_50year_mtn_2%"
##workingDir_LM <- "LM_Projection_50year_mtn_1%"
##workingDir_LM <- "LM_Projection_50year_wlg"
workingDir_IBM <- "IBM_Projection_50year_mtn_3%"
##workingDir_IBM <- "IBM_Projection_50year_mtn_2%"
##workingDir_IBM <- "IBM_Projection_50year_mtn_1%"
##workingDir_IBM <- "IBM_Projection_50year_wlg"

setwd(workingDir_LM)
##setwd(workingDir_IBM)

allScenarioFiles <- list.files(pattern="*.csv")
for (i in 1:length(allScenarioFiles)){
  assign(allScenarioFiles[i], 
         read.csv(allScenarioFiles[i], header=TRUE)
  )}

########## IMPORTANT !!!!!! #########
########## IMPORTANT !!!!!! #########
########## IMPORTANT !!!!!! #########
## Select simulation objects from list. Uncommented depending on whether IBM or LM:
finalPopObjects <- c("Results_LM_Nfinal_MTN_3%.csv","Results_LM_Nfinal_MTN_2%.csv","Results_LM_Nfinal_MTN_1%.csv","Results_LM_Nfinal_WLG.csv","Results_IBM_Nfinal_MTN_3%.csv","Results_IBM_Nfinal_MTN_2%.csv","Results_IBM_Nfinal_MTN_2%.csv","Results_LM_Nfinal_WLG.csv")

startingPopObjects <- c("N40_mtn0.99_IBM.csv", "N40_mtn0.85_IBM.csv", "N40_mtn0.65_IBM.csv")

## lambda is the finite rate of increase of a population over one time step. r is the intrinsinc rate of growth. negative r values indicate a population in decline. lambda < 1 indicates a decline. the relationship between lambda and r : lambda = Nt+1  / Nt, r = ln(lambda), lambda = e^r
growthRates <- data.frame(scenario = as.factor(LETTERS[1:9]), 
                           growth_rate = NA, lambda = NA)
index<-0
for(j in 1:9){
  index <- index+1
  finalPop <- get(finalPopObjects[3]) ## here, j indicates which csv file to read
  startPop <- get(startingPopObjects[3])
  logLambda <- mean((1/10)*log((finalPop[,j])/(startPop[,j]))) ## 50 years for the census time period, loglambda = 1/timeperiod*log(Ntfinal)/Nt0
  lambda <- exp(logLambda)
  growthRates[index,2:3] <- c(logLambda, lambda)
}

growthRates
