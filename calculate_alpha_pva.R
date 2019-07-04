## Select the correct folder for either WLG or MTN data
##workingDir <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_IBM_50year_mtn_0.99/")
##workingDir <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_IBM_50year_mtn_0.85/")
##workingDir <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_IBM_50year_mtn_0.65/")
##workingDir <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_IBM_50year_wlg_0.42/")
##workingDir <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_LM_50year_mtn_3%/")
##workingDir <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_LM_50year_mtn_2%/")
##workingDir <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_LM_50year_mtn_1%/")
##workingDir <- ("/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/pva_LM_50year_wlg/")

setwd(workingDir)
allScenarioFiles <- list.files(pattern="*.csv")
for (i in 1:length(allScenarioFiles)){
  assign(allScenarioFiles[i], 
         read.csv(paste(workingDir, allScenarioFiles[i], sep=''), header=TRUE)
  )}

########## IMPORTANT !!!!!! #########
########## IMPORTANT !!!!!! #########
########## IMPORTANT !!!!!! #########
## Select simulation objects from list. Uncommented depending on whether IBM or LM:
stochObjects <- c("IBM_Scenario1.csv","IBM_Scenario2.csv","IBM_Scenario3.csv","IBM_Scenario4.csv","IBM_Scenario5.csv","IBM_Scenario6.csv","IBM_Scenario7.csv","IBM_Scenario8.csv","IBM_Scenario9.csv")
##stochObjects <- c("LM_Scenario1.csv","LM_Scenario2.csv","LM_Scenario3.csv","LM_Scenario4.csv","LM_Scenario5.csv","LM_Scenario6.csv","LM_Scenario7.csv","LM_Scenario8.csv","LM_Scenario9.csv")

Nfinal <- data.frame(matrix(ncol=9, nrow=1000))
colnames(Nfinal) <- LETTERS[1:9]

index <- 0
for(j in 1:length(stochObjects)){
  index<-index+1
  resX <- get(stochObjects[j])
  Nfinal[,j] <- as.numeric(resX[480,])
}

##write.csv(Nfinal, file=paste0(workingDir_Results,"N40_mtn0.99_IBM.csv"), row.names=F)
##write.csv(Nfinal, file=paste0(workingDir_Results,"N40_mtn0.85_IBM.csv"), row.names=F)
##write.csv(Nfinal, file=paste0(workingDir_Results,"N40_mtn0.65_IBM.csv"), row.names=F)
##write.csv(Nfinal, file=paste0(workingDir_Results,"Nfinal_mtn0.85_IBM.csv"), row.names=F)
##write.csv(Nfinal, file=paste0(workingDir_Results,"Nfinal_mtn0.65_IBM.csv"), row.names=F)
##write.csv(Nfinal, file=paste0(workingDir_Results,"Nfinal_wlg0.42_IBM.csv"), row.names=F)
##write.csv(Nfinal, file=paste0(workingDir_Results,"Nfinal_mtn3%_LM.csv"), row.names=F)
##write.csv(Nfinal, file=paste0(workingDir_Results,"Nfinal_mtn2%_LM.csv"), row.names=F)
##write.csv(Nfinal, file=paste0(workingDir_Results,"Nfinal_mtn1%_LM.csv"), row.names=F)
##write.csv(Nfinal, file=paste0(workingDir_Results,"Nfinal_wlg_LM.csv"), row.names=F)

#######################################################
############### DETERMINE CORRECT alpha ###############
#######################################################

workingDir <- "~/Users/neethaiyer/Desktop/PVA_Tshiaberimu_R/"
setwd(workingDir)
source("1. Function Definitions.R") 

#### Rerun deterministic projection with No=30, t= 50 years
## Reintroduction Scenarios:
ReintroScenario <- read.csv(paste0(workingDir, "ReintroductionScenarios_LM.csv")) ## csv file with Reintroduction Scenarios for LM

mat <- as.matrix(read.csv(paste0(workingDir, "LeslieMatrix_MTN_3%.csv"))) ## csv file with appropriate Leslie Matrix (needs to be converted to matrix object!)

dat <- read.csv(paste0(workingDir, "Gorilla_LifeTables.csv"))
dat$fertilityrate_2percent <- dat[,3]*.786 
## fertility rates multiplied by factor less than 1 to get eigen values of 1.01 which corresponds to a 1% growth rate
dat$fertilityrate_1percent <- dat[,3]*.643 
## fertility rates multiplied by factor less than 1 to get eigen values of 1.02 which corresponds to a 2% growth rate

nyears <- 50 ## Projection Period
nruns <- 1000 ## Number of simulations to run

datX <- dat[,1:3] ## Subset appropriate life history columns: dat[,c(1,4:5)] for WLG, dat[,1:3] for MTN
## NOTE: this subsetting is needed because columns for dat are specified in FUNCTIONS 8 and 9
weaningAge <- 3.5 ## 4.5 for WLG, 3.5 for MTN
adultAge <- 8 ## 10 for WLG, 8 for MTN

temp <- pop_projection(tfinal=nyears, LM=mat, No=ReintroScenario[,11]) ## run the LM projection
Nfinal <- data.frame(ReintroScenario$age, temp[,51])
colnames(Nfinal) <- c("age","N")
sum(Nfinal$N) ## check if Nfinal is about the same as Nfinal_expected for the required lambda value
Nfinal_expected <- 1.032^50*30 ## Nfinal_expected = lambda^t*No
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
alpha <- 0.35 ## set alpha value

res <- matrix(0, nrow=trunc(nyears/timeunit)+1, ncol=nruns)
for(j in 1:length(initalConditions)){
  for(i in 1:nruns){
    print(i)
    abmDataLog <- simTshia(ages0 = initalConditions[[j]][,1], status0 = initalConditions[[j]][,2], time0 = initalConditions[[j]][,3], nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F)
    nindiv <- tapply(abmDataLog$status,abmDataLog$timestep, function(v) length(v)+rbinom(1, sum(v=="L"), .5))##we're adding the unweaned females
    res[1:length(nindiv),i] <- nindiv
  }
  ##write.csv(res, file=paste0("pva_projection_IBM/IBM_Scenario", j,".csv"), row.names=F)
}

finalPop <- as.numeric(res[nrow(res),])
startPop <- as.numeric(res[480,])
##startPop <- rep(100,nruns)

logLambda <- mean((1/10)*log(finalPop/startPop)) ## nyears for the census time period, loglambda = 1/timeperiod*log(Ntfinal)/Nt0
lambda <- exp(logLambda)
growthRates[17,1:3] <- c(alpha, round(logLambda, digits=3), round(lambda, digits=3))
growthRates <- growthRates[order(-growthRates$alpha_value),] 
growthRates

write.csv(growthRates, file=paste0("pva_projection_IBM","growthRates.csv"), row.names=F)
## change row to 15 for next alpha

