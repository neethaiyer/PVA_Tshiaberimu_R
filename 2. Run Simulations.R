## Set the working directory:
workingDir <- "/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/"
setwd(workingDir)
workingDir_Results <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_extn_results/")
##workingDir <- "~/Documents/git repositories/PVA_Tshiaberimu_R/"

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

dat <- read.csv(paste0(workingDir, "Gorilla_LifeTables.csv"))
dat$fertilityrate_2percent <- dat[,3]*.786 
## fertility rates multiplied by factor less than 1 to get eigen values of 1.01 which corresponds to a 1% growth rate
dat$fertilityrate_1percent <- dat[,3]*.643 
## fertility rates multiplied by factor less than 1 to get eigen values of 1.02 which corresponds to a 2% growth rate

####################################################################################
################# OPTIONAL: CREATE CSV FILES with LESLIE MATRICES ##################
####################################################################################

leslieMatrix(lifetable=dat[,1:3], filename=paste0(workingDir,"LeslieMatrix_MTN_3%.csv"))
leslieMatrix(lifetable=dat[,c(1:2, 6)], filename=paste0(workingDir,"LeslieMatrix_MTN_2%.csv"))
leslieMatrix(lifetable=dat[,c(1:2, 7)], filename=paste0(workingDir,"LeslieMatrix_MTN_1%.csv"))
leslieMatrix(lifetable=dat[,c(1, 4:5)], filename=paste0(workingDir,"LeslieMatrix_WLG.csv"))

###############################################################################
############## SET THE INITIAL CONDITIONS OF THE LM & IBM MODELS ##############
###############################################################################

## Reintroduction Scenarios:
ReintroScenario <- read.csv(paste0(workingDir, "ReintroductionScenarios_LM.csv")) ## csv file with Reintroduction Scenarios for LM
ReintroScenario_IBM <- read.csv(paste0(workingDir, "ReintroductionScenarios_IBM.csv")) ## csv file with Reintroduction Scenarios for IBM

## Time parameters:
mat <- as.matrix(read.csv(paste0(workingDir, "LeslieMatrix_WLG.csv"))) ## csv file with appropriate Leslie Matrix (needs to be converted to matrix object!)
nyears <- 50 ## Projection Period
nruns <- 1000 ## Number of simulations to run
timeunit <- 1/12 ## timestep for IBM

## Initial demographic parameters: survivorship, fertility, and weaning age
datX <- dat[,1:3] ## Subset appropriate life history columns: dat[,c(1,4:5)] for WLG, dat[,1:3] for MTN
## NOTE: this subsetting is needed because columns for dat are specified in FUNCTIONS 8 and 9
weaningAge <- 3.5 ## 4.5 for WLG, 3.5 for MTN
adultAge <- 8 ## 10 for WLG, 8 for MTN
alpha <- 0.99 ## function of the fertility rate, 0.99, 0.85 and 0.65 for MTN gorillas with 3%, 2%, and 1% growth rates and 0.42 for WLG growth rates. 

## Depending on the adult female age and weaning age, create a list with the starting conditions for each scenario of the IBM
initalConditions <- convertToList(scenario = ReintroScenario_IBM, adultAge=adultAge, weaningAge=weaningAge)

#####################################################################################
##################### PART 1: Leslie Matrix model (Simple PVA) ######################
#####################################################################################

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
    write.csv(det, file=paste0(workingDir, "pva_projection_LM/LM_Det_Scenario", i-1, ".csv"), row.names=F)
  }

## Second, use the stochastic function:
temp <- matrix(0, nrow=nyears+1, ncol=nruns)
## create an empty matrix that will save the number of individuals for each year of the projection for each run of the LM projection 

for(j in 1:length(ReintroScenario)){
  for(i in 1:nruns) {
    temp[1:(nyears+1),i] <- apply(stoch_projection(tfinal=nyears, LM=mat, No=ReintroScenario[,j+1]),2,sum)
  }
  write.csv(temp, file=paste0(workingDir,"pva_projection_LM/LM_Scenario", j,".csv"), row.names=F)
}

###############################################################################
################ CALCULATE LIKELIHOOD OF EXTINCTION & Pop = 50 ################
###############################################################################
### Now that the csv files have been written, we may not want to re-run the code for as long in the future, so we can just read the generated files and plot the data directly. Make sure you're reading the csv files from the correct folder. 

## Select the correct folder for either WLG or MTN data
workingDir_LM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_LM_50year_mtn_3%/")
##workingDir_LM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_LM_50year_mtn_2%/")
##workingDir_LM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_LM_50year_mtn_1%/")
##workingDir_LM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_LM_50year_wlg/")

setwd(workingDir_LM)
allScenarioFiles <- list.files(pattern="*.csv")

## Calculate the probability that a simulation results in extinction at the end of 50 years. 
## Calculate the probability that the population reaches at least 50 individuals within 50 years.

prob_50years <- data.frame(scenario = as.factor(LETTERS[1:9]), 
                           prob_50 = NA, 
                           prob_Extn = NA)
index <- 0
for(i in 1:length(allScenarioFiles)){
  index <- index+1
  tempx <- as.matrix(read.csv(paste0(workingDir_LM, allScenarioFiles[[i]])))
  ext <- tempx[nrow(tempx),]==0
  probExt <- mean(ext)
  probNe_50 <- mean(tempx[nrow(tempx),]>=50)
  prob_50years[index,2:3] <- c(probNe_50, probExt)*100	
}

## Make sure you use the correct LM for WLG or MTN when you run the code
##dir.create(paste0(workingDir,"pva_extn_results"))##run this line if pva_extn_results does not exist and needs to be created
## MTN projection 3% growth
write.csv(prob_50years, file=paste0(workingDir,"pva_extn_results/extn_lm_MTN_3%.csv"), row.names=F)
## MTN projection 2% growth
##write.csv(prob_50years, file=paste0(workingDir,"pva_extn_results/extn_lm_MTN_2%.csv"), row.names=F)
## MTN projection 1% growth
##write.csv(prob_50years, file=paste0(workingDir,"pva_extn_results/extn_lm_MTN_1%.csv"), row.names=F)
## WLG projection
##write.csv(prob_50years, file=paste0(workingDir,"pva_extn_results/extn_lm_WLG.csv"), row.names=F)

#####################################################################################
################ PART 2: Individual Based Model (IBM) (Complex PVA) #################
#####################################################################################

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
  write.csv(res, file=paste0(workingDir,"pva_projection_IBM/IBM_Scenario", j,".csv"), row.names=F)
}

###############################################################################
################ CALCULATE RISK OF EXTINCTION & Prob Pop = 50 #################
###############################################################################

### Now that the csv files have been written, we may not want to re-run the code for as long in the future, so we can just read the generated files and plot the data directly. Make sure you're reading the csv files from the correct folder. 

## Select the correct folder for either WLG or MTN data
workingDir_IBM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_IBM_50year_mtn_0.99/")
##workingDir_IBM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_IBM_50year_mtn_0.85/")
##workingDir_IBM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_IBM_50year_mtn_0.65/")
##workingDir_IBM <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_IBM_50year_wlg_0.42/")

setwd(workingDir_IBM)
allScenarioFiles <- list.files(pattern="*.csv")

## Calculate the probability that a simulation results in extinction at the end of 50 years. 
## Calculate the probability that the population reaches at least 50 individuals within 50 years.

prob_50years <- data.frame(scenario = as.factor(LETTERS[1:length(allScenarioFiles)]), 
                           probNe_50 = NA, 
                           prob_Extn = NA)

index <- 0
for(i in 1:length(allScenarioFiles)){
  index <- index+1
  resx <- as.matrix(read.csv(paste0(workingDir_IBM, allScenarioFiles[[i]])))
  probNe_50 <- mean(resx[nrow(resx),]>=50)
  prob_Extn <- mean(resx[nrow(resx),]==0)
  prob_50years[index,2:3] <- c(probNe_50, prob_Extn)*100
}

## Make sure you use the correct LM for WLG or MTN when you run the code
#dir.create(paste0(workingDir,"pva_extn_results"))##run this line if pva_extn_results does not exist and needs to be created
## MTN projection 3% growth
write.csv(prob_50years, file=paste0(workingDir,"pva_extn_results/extn_ibm_MTN_3%.csv"), row.names=F)
#write.csv(prob_50years, file=paste0(workingDir,"pva_extn_results/extn_ibm_MTN_2%.csv"), row.names=F)
#write.csv(prob_50years, file=paste0(workingDir,"pva_extn_results/extn_ibm_MTN_1%.csv"), row.names=F)
#write.csv(prob_50years, file=paste0(workingDir,"pva_extn_results/extn_ibm_WLG.csv"), row.names=F)

###############################################################################
################# CREATE CSV FILES FOR FINAL POPULATION SIZES #################
###############################################################################

## Select the correct folder for either WLG or MTN data
##workingDir <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_IBM_50year_mtn_0.99")
##workingDir <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_IBM_50year_mtn_0.85/")
##workingDir <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_IBM_50year_mtn_0.65/")
##workingDir <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_IBM_50year_wlg_0.42/")
workingDir <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_LM_50year_mtn_3%/")
##workingDir <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_LM_50year_mtn_2%/")
##workingDir <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_LM_50year_mtn_1%/")
##workingDir <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_LM_50year_wlg/")

setwd(workingDir)
allScenarioFiles <- list.files(pattern="*.csv")
for (i in 1:length(allScenarioFiles)){
  assign(allScenarioFiles[i], 
         read.csv(paste(workingDir, allScenarioFiles[i], sep=''), header=TRUE)
  )}

## Select simulation objects from list. Uncommented depending on whether IBM or LM:
stochObjects <- c("IBM_Scenario1.csv","IBM_Scenario2.csv","IBM_Scenario3.csv","IBM_Scenario4.csv","IBM_Scenario5.csv","IBM_Scenario6.csv","IBM_Scenario7.csv","IBM_Scenario8.csv","IBM_Scenario9.csv")
##stochObjects <- c("LM_Scenario1.csv","LM_Scenario2.csv","LM_Scenario3.csv","LM_Scenario4.csv","LM_Scenario5.csv","LM_Scenario6.csv","LM_Scenario7.csv","LM_Scenario8.csv","LM_Scenario9.csv")

Nfinal <- data.frame(matrix(ncol=9, nrow=1000))
colnames(Nfinal) <- LETTERS[1:9]

index <- 0
for(j in 1:length(stochObjects)){
  index<-index+1
  resX <- get(stochObjects[j])
  Nfinal[,j] <- as.numeric(resX[nrow(resX),])
}

write.csv(Nfinal, file=paste0(workingDir_Results,"Nfinal_mtn0.99_IBM.csv"), row.names=F)
##write.csv(Nfinal, file=paste0(workingDir_Results,"Nfinal_mtn0.85_IBM.csv"), row.names=F)
##write.csv(Nfinal, file=paste0(workingDir_Results,"Nfinal_mtn0.65_IBM.csv"), row.names=F)
##write.csv(Nfinal, file=paste0(workingDir_Results,"Nfinal_wlg0.42_IBM.csv"), row.names=F)
##write.csv(Nfinal, file=paste0(workingDir_Results,"Nfinal_mtn3%_LM.csv"), row.names=F)
##write.csv(Nfinal, file=paste0(workingDir_Results,"Nfinal_mtn2%_LM.csv"), row.names=F)
##write.csv(Nfinal, file=paste0(workingDir_Results,"Nfinal_mtn1%_LM.csv"), row.names=F)
##write.csv(Nfinal, file=paste0(workingDir_Results,"Nfinal_wlg_LM.csv"), row.names=F)

###############################################################################
################### CALCULATE GROWTH RATE R FOR EACH MODEL ####################
###############################################################################
setwd(workingDir_Results)

allScenarioFiles <- list.files(pattern="*.csv")
for (i in 1:length(allScenarioFiles)){
  assign(allScenarioFiles[i], 
         read.csv(paste(workingDir_Results, allScenarioFiles[i], sep=''), header=TRUE)
  )}

## Select simulation objects from list. Uncommented depending on whether IBM or LM:
finalPopObjects <- c("Nfinal_mtn0.99_IBM.csv","Nfinal_mtn0.85_IBM.csv","Nfinal_mtn0.65_IBM.csv","Nfinal_wlg0.42_IBM.csv","Nfinal_mtn3%_LM.csv","Nfinal_mtn2%_LM.csv","Nfinal_mtn1%_LM.csv","Nfinal_wlg_LM.csv")

## lambda is the finite rate of increase of a population over one time step. r is the intrinsinc rate of growth. negative r values indicate a population in decline. lambda < 1 indicates a decline. the relationship between lambda and r : lambda = Nt+1  / Nt, r = ln(lambda), lambda = e^r
growthRates <- data.frame(scenario = as.factor(LETTERS[1:9]), 
                           growth_rate = NA, lambda = NA)
index<-0
for(j in 1:9){
  index <- index+1
  finalPop <- get(finalPopObjects[1]) ## here, j indicates which csv file to read
  logLambda <- (1/50)*log(index/mean(finalPop[,j])) ## 50 years for the census time period, loglambda = 1/timeperiod*log(Ntfinal)/Nt0
  lambda <- exp(logLambda)
  growthRates[index,2:3] <- c(logLambda, lambda)
}

growthRates
