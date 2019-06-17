## Set the working directory:
workingDir <- "~/Box Sync/PVA_Paper/PVA_Tshiaberimu_R/"
##workingDir <- "~/Documents/git repositories/PVA_Tshiaberimu_R/"

source("1. Function Definition.R")

#####################################################################################
##################### PART 1: Leslie Matrix model (Simple PVA) ######################
#####################################################################################

## Projections were conducted using 4 possible growth rates: 
## A. Mountain gorillas with fertility rates that correspond to 3% growth rate
## B. Mountain gorillas with fertility rates that correspond to 2% growth rate
## C. Mountain gorillas with fertility rates that correspond to 1% growth rate
## D. Western lowland gorillas with constant fertility rates (due to limited data)
## Age of first reproduction = 8 years (for MTN) and 10 years (for WLG)

###############################################################################
#################### SET THE INITIAL CONDITIONS OF THE LM #####################
###############################################################################

## Reintroduction Scenarios for LM
ReintroScenario <- read.csv(paste0(workingDir, "ReintroductionScenarios_LM.csv"))

nyears <- 50 ## Projection Period
nruns <- 1000 ## Number of simulations to run
## Leslie Matrices
matMTN3 <- read.csv(paste0(workingDir, "LeslieMatrix_MTN_3%.csv"))
matMTN2 <- read.csv(paste0(workingDir, "LeslieMatrix_MTN_2%.csv"))
matMTN1 <- read.csv(paste0(workingDir, "LeslieMatrix_MTN_1%.csv"))
matWLG <- read.csv(paste0(workingDir, "LeslieMatrix_WLG.csv"))
mat <- matWLG ## Select the appropriate matrix

###############################################################################
############################# RUN THE LM FUNCTIONS ############################
###############################################################################

## Apply LM projection functions using each reintrodcution scenario
## First, use the deterministic function:
projectPop <- for(i in 2:ncol(ReintroScenario)){
  No <- ReintroScenario[,i] ## Get the reintroduction scenario
  N <- pop_projection(tfinal=nyears, LM=mat, No=No) ## Apply pop_projection function to No for this scenario
  scenario <- strsplit(colnames(ReintroScenario)[i], "_")[[1]][2] ## Get the last element of the column name for each reintroducion scenario
  assign(paste0("N_projected_det", scenario), apply(N,2,sum))  ## Apply the deterministic projection to all scenarios. The assign function takes a variable name as a character string and assigns a value to it. In this case, the values are N at each time step of the projection
}

## Second, use the stochastic function: 
tempA <- tempB <- tempC <- tempD <- tempE <- tempF <- tempG <- tempH <- tempI <- matrix(0, nrow=nyears+1, ncol=nruns)
## create empty matrices that will save the number of individuals for each year of the projection for each run of the LM projection 

## Run each scenario 1000 times using the stochastic projection
for(i in 1:nruns) {
  tempA[1:(nyears+1),i] <- apply(stoch_projection(tfinal=nyears, LM=mat, No=ReintroScenario$No_A),2,sum)
  tempB[1:(nyears+1),i] <- apply(stoch_projection(tfinal=nyears, LM=mat, No=ReintroScenario$No_B),2,sum)
  tempC[1:(nyears+1),i] <- apply(stoch_projection(tfinal=nyears, LM=mat, No=ReintroScenario$No_C),2,sum)
  tempD[1:(nyears+1),i] <- apply(stoch_projection(tfinal=nyears, LM=mat, No=ReintroScenario$No_D),2,sum)
  tempE[1:(nyears+1),i] <- apply(stoch_projection(tfinal=nyears, LM=mat, No=ReintroScenario$No_E),2,sum)
  tempF[1:(nyears+1),i] <- apply(stoch_projection(tfinal=nyears, LM=mat, No=ReintroScenario$No_F),2,sum)
  tempG[1:(nyears+1),i] <- apply(stoch_projection(tfinal=nyears, LM=mat, No=ReintroScenario$No_G),2,sum)
  tempH[1:(nyears+1),i] <- apply(stoch_projection(tfinal=nyears, LM=mat, No=ReintroScenario$No_H),2,sum)
  tempI[1:(nyears+1),i] <- apply(stoch_projection(tfinal=nyears, LM=mat, No=ReintroScenario$No_I),2,sum)
}

## Calculate the probability that a simulation results in extinction at the end of 50 years. 
## Calculate the probability that the population reaches at least 50 (or 40, 100, 150) individuals within 50 years?
prob_50years <- data.frame(scenario = as.factor(LETTERS[1:9]), 
                           prob_150 =  NA, 
                           prob_100 = NA, 
                           prob_50 = NA, 
                           prob_40 = NA, 
                           prob_Extn = NA)
index <- 0
for(i in c("tempA", "tempB", "tempC", "tempD", "tempE", "tempF", "tempG", "tempH", "tempI")){
  index <- index+1
  tempx <- get(i)
  ext <- tempx[nrow(tempx),]==0
  probExt <- mean(ext)
  probNe_150 <- mean(tempx[nrow(tempx),]>=150)
  probNe_100 <- mean(tempx[nrow(tempx),]>=100)
  probNe_50 <- mean(tempx[nrow(tempx),]>=50)
  probNe_40 <- mean(tempx[nrow(tempx),]>=40)
  prob_50years[index,2:6] <- c(probNe_150, probNe_100, probNe_50, probNe_40, probExt)*100	
}

## Make sure you use the correct LM for WLG or MTN when you run the code
#dir.create(paste0(workingDir,"pva_lambda_extn"))##run this line if pva_lambda_extn does not exist and needs to be created
write.csv(prob_50years, file=paste0(workingDir,"pva_lambda_extn/extn_lm_WLG.csv"), row.names=F)
## MTN projection 3% growth
write.csv(prob_50years, file=paste0(workingDir,"pva_lambda_extn/extn_lm_MTN.csv"), row.names=F)
## MTN projection 2% growth
write.csv(prob_50years, file=paste0(workingDir,"pva_lambda_extn/extn_lm_MTN_1percent.csv"), row.names=F)
## MTN projection 1% growth
write.csv(prob_50years, file=paste0(workingDir,"pva_lambda_extn/extn_lm_MTN_2percent.csv"), row.names=F)

#####################################################################################
################ PART 2: Individual Based Model (IBM) (Complex PVA) #################
#####################################################################################

###############################################################################
#################### SET THE INITIAL CONDITIONS OF THE IBM ####################
###############################################################################

## Reintroduction Scenarios for IBM
ReintroScenario_IBM <- read.csv(paste0(workingDir, "ReintroductionScenarios_IBM.csv"))

## Depending on the adult female age and weaning age, create a list with the starting conditions for each scenario of the IBM
initalConditions <- convertToList(scenario = ReintroScenario_IBM, adultAge=10, weaningAge=4.5)
nyears <- 50 ## Projection Period
timeunit <- 1/12 ## timestep
nruns <- 1000 ## Number of simulations to run
alpha <- 0.42 ## function of the fertility rate, 0.99 for MTN 3% growth rate, 0.42 for WLG, 

###############################################################################
############################# RUN THE IBM FUNCTION ############################
###############################################################################

res <- matrix(0, nrow=trunc(nyears/timeunit)+1, ncol=nruns)
for(j in 1:length(initalConditions)){
  for(i in 1:nruns){
    print(i)
    abmDataLog <- simTshia(ages0 = initalConditions[[j]][,1], status0 = initalConditions[[j]][,2], time0 = initalConditions[[j]][,3], nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F)
    nindiv <- tapply(abmDataLog$status,abmDataLog$timestep, function(v) length(v)+rbinom(1, sum(v=="L"), .5))##we're adding the unweaned females
    res[1:length(nindiv),i] <- nindiv
  }
  write.csv(res, file=paste0(workingDir,"pva_IBM_50year/Scenario", j,".csv"), row.names=F)
}
