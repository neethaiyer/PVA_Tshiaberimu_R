## Set the working directory:
workingDir <- "~/Box Sync/PVA_Paper/PVA_Tshiaberimu_R/"
##workingDir <- "~/Documents/git repositories/PVA_Tshiaberimu_R/"

source("1. Function Definition.R")

#####################################################################################
############################## LOAD ALL CSV FILES ###################################
#####################################################################################

## Mountain gorilla (MTN) life history parameters taken from Bronikowski et al (2016)
## Western lowland gorilla (WLG) life history parameters taken from Breuer et al (2010), Breuer (2008) 
dat <- read.csv(paste0(workingDir, "Gorilla_LifeTables.csv"))
dat$fertilityrate_2percent <- dat[,3]*.786 
## fertility rates multiplied by factor less than 1 to get eigen values of 1.01 which corresponds to a 1% growth rate
dat$fertilityrate_1percent <- dat[,3]*.643 
## fertility rates multiplied by factor less than 1 to get eigen values of 1.02 which corresponds to a 2% growth rate

## Reintroduction Scenarios for LM
ReintroScenario <- read.csv(paste0(workingDir, "ReintroductionScenarios_LM.csv"))

#####################################################################################
##################### PART 1: Leslie Matrix model (Simple PVA) ######################
#####################################################################################

## Projections were conducted using 4 possible growth rates: 
## A. Mountain gorillas with fertility rates that correspond to 3% growth rate
## B. Mountain gorillas with fertility rates that correspond to 2% growth rate
## C. Mountain gorillas with fertility rates that correspond to 1% growth rate
## D. Western lowland gorillas with constant fertility rates (due to limited data)
## Age of first reproduction = 8 years (for MTN) and 10 years (for WLG)

########################################################################################
#################### STEP 1: CREATE LESLIE MATRICES FROM LIFE TABLES ###################
########################################################################################

## A. Create the Leslie Matrix for MTN with fertility rates that correspond to 3% growth rate
mat <- matrix(0, nrow=nrow(dat), ncol=nrow(dat)) ## create an empty square matrix
mat[1,] <- dat[,3] ## first row in matrix assigned the fertility rates from the life table
mat2 <- matrix(0,ncol=ncol(mat)-1, nrow=ncol(mat)-1) 
diag(mat2) <- 1-dat[-nrow(dat),2] ## survival rates are assigned to just under the diagonal of a LM
mat[2:nrow(mat), 1:(ncol(mat)-1)] <- mat2
write.csv(mat, file=paste0(workingDir,"LeslieMatrix_MTN_3%.csv"), row.names=F)

## B. Create the Leslie Matrix for MTN with fertility rates that correspond to 2% growth rate
mat <- matrix(0, nrow=nrow(dat), ncol=nrow(dat)) ## create an empty square matrix
mat[1,] <- dat[,6] ## first row in matrix assigned the fertility rates from the life table
mat2 <- matrix(0,ncol=ncol(mat)-1, nrow=ncol(mat)-1) 
diag(mat2) <- 1-dat[-nrow(dat),2] ## survival rates are assigned to just under the diagonal of a LM
mat[2:nrow(mat), 1:(ncol(mat)-1)] <- mat2
write.csv(mat, file=paste0(workingDir,"LeslieMatrix_MTN_2%.csv"), row.names=F)

## C. Create the Leslie Matrix for MTN with fertility rates that correspond to 1% growth rate
mat <- matrix(0, nrow=nrow(dat), ncol=nrow(dat)) ## create an empty square matrix
mat[1,] <- dat[,7] ## first row in matrix assigned the fertility rates from the life table
mat2 <- matrix(0,ncol=ncol(mat)-1, nrow=ncol(mat)-1) 
diag(mat2) <- 1-dat[-nrow(dat),2] ## survival rates are assigned to just under the diagonal of a LM
mat[2:nrow(mat), 1:(ncol(mat)-1)] <- mat2
write.csv(mat, file=paste0(workingDir,"LeslieMatrix_MTN_1%.csv"), row.names=F)

## D. Create the Leslie Matrix for WLG with constant fertility rates
mat <- matrix(0, nrow=nrow(dat), ncol=nrow(dat)) ## create an empty square matrix
mat[1,] <- dat[,5] ## first row in matrix assigned the fertility rates from the life table
mat2 <- matrix(0,ncol=ncol(mat)-1, nrow=ncol(mat)-1) 
diag(mat2) <- 1-dat[-nrow(dat),4] ## survival rates are assigned to just under the diagonal of a LM
mat[2:nrow(mat), 1:(ncol(mat)-1)] <- mat2
write.csv(mat, file=paste0(workingDir,"LeslieMatrix_WLG.csv"), row.names=F)

########################################################################################
########################## STEP 2: CREATE DEMOGRAPHIC PYRAMIDS #########################
########################################################################################

## A. Demographic pyramid for MTN (same for 3%, 2%, and 1%)
n <- rep(1, nrow(dat))
n[1] <- 1
for (i in 2:length(n)){
  n[i] <- prod(1-dat[1:(i-1),2])
} ## Make sure sum equals 1 to generate pyramid
n_mtn <- n/(sum(n))
## for the cumulative survival curve:
n_mtnCS <- n/n[1]

## D. Demographic pyramid for WLG
n <- rep(1, nrow(dat))
n[1] <- 1
for (i in 2:length(n)){
  n[i] <- prod(1-dat[1:(i-1),4])
} ## Make sure sum equals 1 to generate pyramid
n_wlg <- n/(sum(n))
## for the cumulative survival curve:
n_wlgCS <- n/n[1]

########################################################################################
######################### STEP 2: RUN LESLIE MATRIX SIMULATIONS ########################
########################################################################################

##############################################
########### SET INITIAL CONDITIONS ###########
##############################################

nyears <- 50 ## projection period
nruns <- 1000 ## number of simulations of the model
## Leslie Matrices
matMTN3 <- read.csv(paste0(workingDir, "LeslieMatrix_MTN_3%.csv"))
matMTN2 <- read.csv(paste0(workingDir, "LeslieMatrix_MTN_2%.csv"))
matMTN1 <- read.csv(paste0(workingDir, "LeslieMatrix_MTN_1%.csv"))
matWLG <- read.csv(paste0(workingDir, "LeslieMatrix_WLG.csv"))

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
