## Set the working directory:
workingDir <- "~/Box Sync/PVA_Paper/PVA_Tshiaberimu_R/"
##workingDir <- "~/Documents/git repositories/PVA_Tshiaberimu_R/"

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

############################################################################################
############################# STEP 3: LESLIE MATRIX PROJECTIONS ############################
############################################################################################

## The function below returns the deterministic LM population size projection, using the projection period (tfinal), age and number of females in the starting population (No) and Leslie matrix (LM). 
#mat <- mat_wlg
pop_projection <- function(tfinal, LM=mat, No=No){
  pop <- N <- No
  for (i in 1:tfinal){
    pop <- LM%*%pop
    N <- cbind(N,pop)
  }
  N
}

## The function below returns the stochastic LM population size projection, using the projection period (tfinal), age and number of females in the starting population (No) and Leslie matrix (LM). 
stoch_projection <- function(tfinal, LM=mat, No=No){
  pop <- N <- No
  for (i in 1:tfinal){
    currentpop <- pop
    for (j in 1:length(No)){
      pop[j] <- sum(rbinom(n=length(No), size=currentpop, prob=LM[j,]))
    }
    N <- cbind(N, pop)
  }
  N
}

#####################################################################################
################ PART 2: Individual Based Model (IBM) (Complex PVA) #################
#####################################################################################

## For this model, we specified 5 categories of individuals: 
## immature (I): either weaned juveniles or unweaned infants, cycling adult female (C), pregnant adult female (P), lactating adult female (L), and females that lose their infants (CD)
## Thus, our IBM will inclde both an age and reproductive status category vector. See STEPS 1-2 below.
## Make sure to choose between either MTN or WLG parameters for the steps below. 

#################################################################################
######################## GORILLA LIFE HISTORY PARAMETERS ########################
#################################################################################

## Initial parameters: survivorship, fertility, and weaning age
datX <- dat[,c(1,4:5)] ## Subset appropriate life history columns: dat[,c(1,4:5)] for WLG, dat[,1:3] for MTN
## ***NOTE: this subsetting is needed because columns for dat are specified in STEP 3 (in 4. and 5.)
weaningAge <- 4.5 ## 4.5 for WLG, 3.5 for MTN
adultAge <- 10 ## 10 for WLG, 8 for MTN
n_pyramid <- n_mtn ## n_wlg or n_mtn
alpha <- 0.42 ## 0.42 for WLG OR 0.99 for MTN

########################################################################################
######################## STEP 1: AGE CATEGORIES FOR INDIVIDUALS ########################
########################################################################################

ages0 <- numeric(0) ## we need a vector that keeps track of the age of individuals
ages0 <- rep(1:length(n_pyramid), rmultinom(1,100,prob=n_pyramid)) ## We start with 100 individuals. rmultinom creates a vector that assigns the number of individuals (in this case 100) to each age, according to the LM demography pyramid
hist(ages0)

## A. Parameters for optimistic model (MTN)
## Set up a variable that keeps track of the time since entry into the current age category
time0 <- numeric(length(ages0))

## Given an individual's age0 at time0, how many years has it been since entering the current age category? 
time0[(ages0-weaningAge)>=0 & (ages0-weaningAge)<weaningAge+1] <- ages0[(ages0-weaningAge)>=0 & (ages0-weaningAge)<weaningAge+1]-weaningAge 
time0[(ages0-adultAge)>=0] <- 0
time0 <- time0[ages0>=weaningAge]
ages0 <- ages0[ages0>=weaningAge]

########################################################################################
######################## STEP 2: BREEDING STATUS OF INDIVIDUALS ########################
########################################################################################

## Adult female status 
status0 <- character(length(ages0))
status0[ages0<adultAge] <- "I" ## Immature individuals are those between 0-8 or 10 years old and non-reproductive
status0[ages0>=adultAge] <- "C" ## We start the simulation with cycling females only. With each timestep in the model, C females may transition to other categories, detailed in Step 3
data.frame(ages0, time0, status0)

########################################################################################
#################### STEP 3: TRANSITION PROBABILITES OF INDIVIDUALS ####################
########################################################################################

## Next, we specifiy the transition probabilities, in relation to time since entry in the "current" category

## Transition probabilities (MTN or WLG)
## t is the time passed in the initial class
timeunit <- 1/12 ## in years. Each time step is 1 month.

## 1. Transition from I to C:
IC <- function(t) ifelse(t<weaningAge+1, 0, 1) ## Probability for an subadult to transition to adult is zero if time as subadult is t<5.5 and 1 if t>=5.5

## 2. Transition from C to P:
alpha <- alpha ## Probability to be pregnant after 12 months being cycling
## relationship between alpha and p: 
## (1-p)^(1/timeunit)=1-alpha
## log(1-p)=log(1-alpha)/(1/timeunit)
## p=1-exp(log(1-alpha)/(1/timeunit))
CP <- function(t, alpha) 1-exp(log(1-alpha)/(1/timeunit)) ## Probability for a cycling adult female to become pregnant, per time step, given the alpha value

## 3. Transition from P to L:
PL <- function(t) ifelse(t<(adultAge+0.5)/12, 0, 1) ## Probability for a pregnant adult female to transition to "lactating" is 1 if she has been pregnant for 8.5 months

## 4. Transition from L to C: (due to infant death)
LCdeathInf <- function(t) 1-exp(log(1-datX[trunc(t+1),2])/(1/timeunit))  ## Probability that a lactating female transitions to cycling due to the loss of her dependent infant
LC <- function(t) ifelse(t<weaningAge,0,1) ## Probability that a lactating female transitions to cycling is 1 after the weaning age of her infant.

## 5. Mortality rate
deathRate <- function(age) 1-exp(log(1-datX[trunc(age+1),2])/(1/timeunit)) ## Death rate, per month, specified by the LM

## Function that specifies the status change of each agent in model : I,C,P,L,CD
statusChange <- function(status, t, alpha){
  if(status =="I"){
    if(rbinom(1,1,IC(t))){
      return("C")
    } else return('I')
  } else if (status =="C"){
    if(rbinom(1,1,CP(t, alpha))){
      return("P")
    } else return('C')
  } else if (status =="P") {
    if(rbinom(1,1,PL(t))){
      return("L")
    } else return('P')
  } else if(!rbinom(1,1,LC(t))){		  ## If her infant doesn't get weaned because it is too young
    if(!rbinom(1,1,LCdeathInf(t))) {	  ## If her infant doesn't die
      return("L")						  ## Then the female is still lactating
    } else return("CD")					  ## If her infant dies, female becomes "CD", i.e. she transitions to cycling after the death of her infant. We distinguish between CD and C because if the baby dies, we do not add it into the population in the next year
  } else  return("C")
}

## Possible statuses: I,C,P,L,CD
## quick test:
## statusChange("L", 0.6, 0.99)
## statusChange("C", 3/12, 0.99)

########################################################################################
#################### STEP 4: SIMULATION FUNCTION FOR STOCHASTIC IBM ####################
########################################################################################

timeunit <- 1/12
alpha <- alpha

simTshia <- function(ages0, status0, time0, nyears=50, timeunit=1/12, alpha=alpha, verbose=T){
  abmDataLog <- data.frame(timestep=0, ages=ages0, time=time0, status=status0, indiv=1:length(ages0), stringsAsFactors = FALSE)
  iter <- max(abmDataLog$indiv)+1 ## iter keeps track of the index of the last row (of abmDataLog) + 1, which is essentially the next available gorilla "ID" in the model. We don't want to repeat ID's. 
  for(i in 1:trunc(nyears/timeunit)){
    if(i%%(1/timeunit)==0 & verbose) print(paste("time =",i*timeunit)) ## %% or modulo find the remainder: every year, when you divide 1 by 12 %% = 0, so this keeps track of time in the simulation. if we select verbose=F, the time is not printed. 
    newAbmData <- abmDataLog[0,]
    abmData <- abmDataLog[abmDataLog$timestep==(i-1),] ## abmData keeps tracks of data from previous timestep
    for(j in 1:nrow(abmData)){ ## j is the focal individual
      currentStatus <- abmData[j,4] ## this is the status of the individual j
      newStatus <- statusChange(currentStatus, abmData[j,3], alpha=alpha) ## Returns "I", "P", "L", "C", or "CD" and newStatus is the status at the next time step
      if(rbinom(1,1,deathRate(abmData[j,2]))==1){ ## Did the individual just die?
        newStatus <- "D"
      }
      if(currentStatus!=newStatus) {
        newtime <- timeunit ## this will depend on the newStatus
        if(newStatus=="L") infantSex <- sample(0:1,1) ## this keeps track of sex, using sex ratio 1:1, where 0=Male, 1=Female
      } else {
        newtime <- abmData[j,3]+timeunit
      }
      newAbmData <- rbind(newAbmData, data.frame(timestep=i, ages=abmData[j,2]+timeunit, time=newtime, status=newStatus, indiv=abmData[j,5], stringsAsFactors = FALSE))
      if(currentStatus=="L" & newStatus=="C"){ ## 50/50 sex ratio
        if( infantSex==1) {
          newAbmData <- rbind(newAbmData, data.frame(timestep=i, ages=weaningAge+timeunit, time=weaningAge+timeunit, status="I", indiv=iter, stringsAsFactors = FALSE))
          iter <- iter+1
        }
      } ## If weaning, add new row to abmData
    }
    if(sum(newAbmData$status=="CD")>0) newAbmData$status[newAbmData$status=="CD"] <- "C"
    if(nrow(newAbmData[newAbmData$status!="D",])==0) break ## if the adult female died, newAbmData should have 1 row with status D, and the juvenile gets weaned, so we keep the baby in the next time step. 
    abmDataLog <- rbind(abmDataLog, newAbmData[newAbmData$status!="D",]) ## if there are individuals that are not dead then add them to the newAbmData
  }
  return(abmDataLog)
}

########################################################################################
######################### STEP 5: SET REINTRODUCTION SCENARIOS #########################
########################################################################################

## Create a list of these scenarios with the age of individuals introduced, and their status at the start of the model
convertToList <- function(scenario, adultAge, weaningAge){
  res <- list()
  for(i in 1:ncol(scenario)){
    ages <- c(na.omit(scenario[,i]))
    ages <- ages[ages>= weaningAge]
    statuses <- ifelse(ages<adultAge, "I", "C")
    times <- numeric(length(ages))
    times <- ages-weaningAge
    times[ages>=adultAge] <- ages[ages>=adultAge]-adultAge
    res[[i]] <- data.frame(ages0 = ages, status0 = statuses, time0 = times)
  }
  return(res)
}
