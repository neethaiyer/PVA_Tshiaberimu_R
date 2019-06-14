## Set the working directory:
workingDir <- "~/Box Sync/PVA_Paper/PVA_Tshiaberimu_R/"
##workingDir <- "~/Documents/git repositories/PVA_Tshiaberimu_R/"

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
datBron <- read.csv(paste0(workingDir, "LifeTable_MTN_Bronikowski.csv"))
dat1 <- datBron
dat1$fertilityrate_1percent <- datBron[,3]*.643 ## fertility rates multiplied by factor less than 1 to get eigen values of 1.01 which corresponds to a 1% growth rate
dat1$fertilityrate_2percent <- datBron[,3]*.786 ## fertility rates multiplied by factor less than 1 to get eigen values of 1.02 which corresponds to a 2% growth rate

## Transform MTN life table into a Leslie matrix (fertility rates are in the top row, survival rates are just under the diagonal)
mat_mtn <- matrix(0, nrow=nrow(dat1), ncol=nrow(dat1)) ## create square matrix with 0s everywhere
mat_mtn[1,] <- dat1[,3] ## first row in matrix assigned the fertility rates from the life table
mat_mtn_2 <- matrix(0,ncol=ncol(mat_mtn)-1, nrow=ncol(mat_mtn)-1) 
diag(mat_mtn_2) <- 1-dat1[-nrow(dat1),2] ## survival rates are assigned to just under the diagonal of a LM
mat_mtn[2:nrow(mat_mtn), 1:(ncol(mat_mtn)-1)] <- mat_mtn_2
head(mat_mtn) ## View and check matrix
#write.csv(mat_mtn, file=paste0(workingDir,"LeslieMatrix_MTN.csv"), row.names=F)
## Calculate the eigenvalue of the matrices
eigenvalues_mtn <- eigen(mat_mtn, only.values=TRUE)
Re(eigenvalues_mtn$values[1]) ## this is the dominant eigenvalue of the MTN LM, i.e. lambda = 1.032567

## Transform MTN life table into a Leslie matrix (fertility rates are in the top row, survival rates are just under the diagonal) : USING GROWTH RATES OF 1%
mat_mtn <- matrix(0, nrow=nrow(dat1), ncol=nrow(dat1)) ## create square matrix with 0s everywhere
mat_mtn[1,] <- dat1[,4] ## first row in matrix assigned the fertility rates from the life table
mat_mtn_2 <- matrix(0,ncol=ncol(mat_mtn)-1, nrow=ncol(mat_mtn)-1) 
diag(mat_mtn_2) <- 1-dat1[-nrow(dat1),2] ## survival rates are assigned to just under the diagonal of a LM
mat_mtn[2:nrow(mat_mtn), 1:(ncol(mat_mtn)-1)] <- mat_mtn_2
head(mat_mtn) ## View and check matrix
mat_mtn_1per <- mat_mtn
#write.csv(mat_mtn, file=paste0(workingDir,"LeslieMatrix_MTN_1percent.csv"), row.names=F)
## Calculate the eigenvalue of the matrices
eigenvalues_mtn <- eigen(mat_mtn, only.values=TRUE)
Re(eigenvalues_mtn$values[1]) ## this is the dominant eigenvalue of the MTN LM, i.e. lambda = 1.010051

## Transform MTN life table into a Leslie matrix (fertility rates are in the top row, survival rates are just under the diagonal) : USING GROWTH RATES OF 2%
mat_mtn <- matrix(0, nrow=nrow(dat1), ncol=nrow(dat1)) ## create square matrix with 0s everywhere
mat_mtn[1,] <- dat1[,5] ## first row in matrix assigned the fertility rates from the life table
mat_mtn_2 <- matrix(0,ncol=ncol(mat_mtn)-1, nrow=ncol(mat_mtn)-1) 
diag(mat_mtn_2) <- 1-dat1[-nrow(dat1),2] ## survival rates are assigned to just under the diagonal of a LM
mat_mtn[2:nrow(mat_mtn), 1:(ncol(mat_mtn)-1)] <- mat_mtn_2
head(mat_mtn) ## View and check matrix
mat_mtn_2per <- mat_mtn
##write.csv(mat_mtn, file=paste0(workingDir,"LeslieMatrix_MTN_2percent.csv"), row.names=F)
## Calculate the eigenvalue of the matrices
eigenvalues_mtn <- eigen(mat_mtn, only.values=TRUE)
Re(eigenvalues_mtn$values[1]) ## this is the dominant eigenvalue of the MTN LM, i.e. lambda = 1.020043

## Transform WLG life table into a Leslie matrix (fertility rates are in the top row, survival rates are just under the diagonal)
mat_wlg <- matrix(0, nrow=nrow(dat), ncol=nrow(dat)) ## create square matrix with 0s everywhere
mat_wlg[1,] <- dat[,3] ## first row in matrix assigned the fertility rates from the life table
mat2 <- matrix(0,ncol=ncol(mat_wlg)-1, nrow=ncol(mat_wlg)-1) 
diag(mat2) <- 1-dat[-nrow(dat),2] ## survival rates are assigned to just under the diagonal of a LM
mat_wlg[2:nrow(mat_wlg), 1:(ncol(mat_wlg)-1)] <- mat2
head(mat_wlg) ## View and check matrix
write.csv(mat_wlg, file=paste0(workingDir,"LeslieMatrix_WLG.csv"), row.names=F)
## Calculate the eigenvalue of the matrices
eigenvalues_wlg <- eigen(mat_wlg, only.values=TRUE)
Re(eigenvalues_wlg$values[1]) ## this is the dominant eigenvalue of the WLG LM, i.e. lambda = 1.020623 vs. updated labmda = 0.9994444

## Let's create a demographic pyramid for WLG
n <- rep(1, nrow(dat))
n[1] <- 1
for (i in 2:length(n)){
  n[i] <- prod(1-dat[1:(i-1),2])
}
## Make sure sum equals 1 to generate pyramid
n <- n/(sum(n))
## for the cumulative survival curve:
n1 <- n/n[1]

## Let's create a demographic pyramid for MTN
n_mtn <- rep(1, nrow(dat1))
n_mtn[1] <- 1
for (i in 2:length(n_mtn)){
  n_mtn[i] <- prod(1-dat1[1:(i-1),2])
}
## Make sure sum equals 1 to generate pyramid
n_mtn <- n_mtn/(sum(n_mtn))
## for the cumulative survival curve:
n_mtn_1 <- n_mtn/n_mtn[1]

## Let's pick our colors for MTN and WLG and include a transparency factor; note that these are only used for the bargraph below:
coral <- rgb(250, 128, 114, alpha=150, maxColorValue = 255) ## for WLG
azure4 <- rgb(131, 139, 139, alpha=150, maxColorValue = 255) ## for MTN

#par(mfrow=c(1,2), oma=c(0,0,0,0), mar=c(5,4,2,1))
## demographic pyramid for MTN
barplot(n_mtn, horiz=T, names.arg=paste0(0:(length(n_mtn)-1), "-", 1:length(n_mtn)), las=1, xlab="relative frequency", col=azure4, cex.axis = 1, cex.names = 0.7, ylab="Age", cex.lab=1, font.lab=2, xlim=c(0,0.07))

## demographic pyramid for WLG
barplot(n, horiz=T, names.arg=paste0(0:(length(n)-1), "-", 1:length(n)), las=1, xlab="relative frequency", col=coral, cex.axis = 1, cex.names = 0.7, ylab="Age", cex.lab=1, font.lab=2, add=TRUE)

## Let's take a look at the annual mortality and cumulative survival curves for MTN and WLG
## 2-pannel survival plot
par(mfrow=c(1,2), oma=c(0,0,0,0), mar=c(5,4,2,1))
plot(dat[,1:2], type="o", bty="l", xlab="Age", ylab="Annual mortality", las=1, bg="coral", cex.lab=0.8, cex.axis=0.8, ylim=c(0,1), cex=0.8, cex.lab=0.8, font.lab=2, pch=21)
lines(dat1[,1:2], type="o", bty="l", xlab="Age", ylab="Annual mortality", las=1, bg="azure4", cex.lab=0.8, cex.axis=0.8, ylim=c(0,1), cex=0.8, cex.lab=0.8, font.lab=2, pch=24)
legend(0, 1, legend=c("Western Gorillas", "Mountain Gorillas"),
       pt.bg=c("coral", "azure4"), lty=c(1,1), cex=0.8, text.font=2, pch=c(21,24))
plot(dat[,1],n1, bty="l", type="o", xlab="Age", ylab="Cumulative survival", las=1, bg="coral", cex.lab=0.8, cex.axis=0.8, ylim=c(0,1), cex=0.8, cex.lab=0.8, font.lab=2, pch=21)
lines(dat1[,1],n_mtn_1, bty="l", type="o", xlab="Age", ylab="Cumulative survival", las=1, bg="azure4", cex.lab=0.8, cex.axis=0.8, ylim=c(0,1), cex=0.8, cex.lab=0.8, font.lab=2, pch=24)
legend(15, 1, legend=c("Western Gorillas", "Mountain Gorillas"),
       pt.bg=c("coral", "azure4"), lty=c(1,1), cex=0.8, text.font=2, pch=c(21,24))

############################################################################################
################################ LM Parameter Specification ################################
############################################################################################

## The asymptotic growth rate lambda was derived from the Leslie matrix (LM). Lambda is the dominant eigenvalue of LM
## At each 1-year time step, we multiplied the Leslie matrix by the current age-specific population estimates to calculate the expected new age-specific population estimates: Nt+1 = LM * Nt
## A stochastic version of this LM model was also created.
## Some assumptions were made when using this LM model: 
## (1) stable environment
## (2) age-specific mortality and fertility rates are constant over time
## (3) population is in the exponential growth phase (no density dependence)
## (4) 1:1 sex ratio at birth
## (5) only females considered and adult females are always reproductively active

## Tshiaberimu Reintroduction Scenarios for simple LM projection
## The dynamics of the Tshiaberimu population under Scenarios A-I specified below was simulated over a timeframe of 50 years
## Scenario A: starting population as current Tshiaberimu status: 1 adult female aged 19 and assuming juvenile is M
## Scenario B: starting population as current Tshiaberimu status: 1 adult female aged 19 and assuming juvenile is F and aged 5
## Scenario C: starting population 1F plus 2 Fs reintroduced / age:7,7,19
## Scenario D: starting population 1F plus 3 Fs reintroduced / age:7,7,8,19
## Scenario E: starting population 1F plus 4 Fs reintroduced / age:7,7,8,9,19
## Scenario F: starting population 1F plus 5 Fs reintroduced / age:7,7,8,9,12,19
## Scenario G: starting population 1F plus 6 Fs reintroduced / age:7,7,8,9,12,12,19
## Scenario H: starting population 1F plus 7 Fs reintroduced / age:7,7,8,9,12,12,17,19
## Scenario I: starting population 1F plus 8 Fs reintroduced / age:7,7,8,9,12,12,17,17,19
  
## Notes:
## 1. Scenarios A & B take into account the uncertainty in sex of the surviving juvenile in the current Tshiaberimu population. Scenarios C-I assume that the current adult female is 19 years and conservatively assume that the juvenile is a male.
## 2. In addition to the deterministic population projections, the mean stochastic population projection and probability of extinction were calculated for 1000 50-year-long simulations.
## 3. Probability of extinction for each scenario was calculated as the total proportion of simulations that resulted in a final population size of 0 females. 

## The function below returns the deterministic LM population size projection, using the projection period (tfinal), age and number of females in the starting population (No) and Leslie matrix (LM). 
mat <- mat_wlg ## make sure you choose either the WLG or MTN LM (mat_wlg or mat_mtn or mat_mtn_1per or mat_mtn_2per respectively)
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

## Now, let's apply the functions using each reintrodcution scenario. Load the appropriate csv file.
ReintroScenario <- read.csv(paste0(workingDir, "ReintroductionScenarios_LM.csv"))

## First, let's use the deterministic function:
nyears <- 50 ## projection period
projectPop <- for(i in 2:ncol(ReintroScenario)){
    No <- ReintroScenario[,i] ## Get the reintroduction scenario
    N <- pop_projection(tfinal=nyears, LM=mat, No=No) ## Apply pop_projection function to No for this scenario
    scenario <- strsplit(colnames(ReintroScenario)[i], "_")[[1]][2] ## Get the last element of the column name for each reintroducion scenario
    assign(paste0("N_projected_det", scenario), apply(N,2,sum))  ## Apply the deterministic projection to all scenarios. The assign function takes a variable name as a character string and assigns a value to it. In this case, the values are N at each time step of the projection
  }

## Second, let's use the stochastic function: 
nruns <- 1000 ## number of simulations of the model
## The temp vectors are empty matrices that will save the number of individuals for each year of the projection for each run of the projection 
tempA <- tempB <- tempC <- tempD <- tempE <- tempF <- tempG <- tempH <- tempI <- matrix(0, nrow=nyears+1, ncol=nruns)

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

  ## stil working here ##
##for(j in 2:ncol(ReintroScenario)){
  No <- ReintroScenario[,j]
  scenario <- strsplit(colnames(ReintroScenario)[j], "_")[[1]][2] ## Get the last element of the column name for each reintroducion scenario
  for(i in 1:nruns){
  assign(paste0("temp", scenario), apply(stoch_projection(tfinal=nyears, LM=mat, No=No),2,sum))
  }
##}
  ## stil working here ##

## What is probability that simulation results in extinction at the end of 50 years? What is the probability that the population reaches at least 50 (or 40, 100, 150) individuals within 50 years?
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
## WLG values
#dir.create(paste0(workingDir,"pva_lambda_extn"))##run this line if pva_lambda_extn does not exist and needs to be created
write.csv(prob_50years, file=paste0(workingDir,"pva_lambda_extn/extn_lm_WLG.csv"), row.names=F)
## MTN values
write.csv(prob_50years, file=paste0(workingDir,"pva_lambda_extn/extn_lm_MTN.csv"), row.names=F)
## MTN values 1%
write.csv(prob_50years, file=paste0(workingDir,"pva_lambda_extn/extn_lm_MTN_1percent.csv"), row.names=F)
## MTN values 2%
write.csv(prob_50years, file=paste0(workingDir,"pva_lambda_extn/extn_lm_MTN_2percent.csv"), row.names=F)

################################################################################################
###################################### PLOTS for LM ############################################
################################################################################################
## You can save the plots using the commented out lines of code below, or directly by saving the plots using RStudio's interface. 
## file_path <- file.path("~.png",paste("PVA",i, ".png", sep = ""))
## png(file_path, width=4, height=4, units='in', res=100)

## Let's plot lambda and extinction probabilities. Use the correct dataframe for either WLG or MTN:
prob_50years_wlg_lm <- read.csv(paste0(workingDir,"pva_lambda_extn/extn_lm_WLG.csv"))
prob_50years_mtn_lm <- read.csv(paste0(workingDir,"pva_lambda_extn/extn_lm_MTN.csv"))
prob_50years_mtn_1per_lm <- read.csv(paste0(workingDir,"pva_lambda_extn/extn_lm_MTN_1percent.csv"))
prob_50years_mtn_2per_lm <- read.csv(paste0(workingDir,"pva_lambda_extn/extn_lm_MTN_2percent.csv"))

################################################################################################
################################ Plots of extinction probability  ##############################
################################################################################################

plot.new()
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,5.1))
plot.window(xlim=c(1,9), ylim=c(0,100))
axis(1, 1:9, LETTERS[1:9])
axis(2)
axis(2, font.lab=2, at=seq(0, 100, by=10), labels=seq(0, 100, by=10))
title(xlab="Reintroduction Scenario", ylab="Probability of Extinction", font.lab=2)
lines(prob_50years_wlg_lm$scenario, prob_50years_wlg_lm$prob_Extn, bg="coral", type="b", pch=21)
lines(prob_50years_mtn_lm$scenario, prob_50years_mtn_lm$prob_Extn, bg="azure4", type="b", pch=24)
lines(prob_50years_mtn_1per_lm$scenario, prob_50years_mtn_1per_lm$prob_Extn, bg="dodgerblue", type="b", pch=22)
lines(prob_50years_mtn_2per_lm$scenario, prob_50years_mtn_2per_lm$prob_Extn, bg="yellow", type="b", pch=23)
legend(1, 100, legend=c("Western Gorillas", "Mountain Gorillas", "1 percent growth", "2 percent growth"),
       pt.bg=c("coral", "azure4","dodgerblue","yellow"), lty=c(1,1,1,1), text.font=2, pch=c(21,24,22,23))

################################################################################################
########################### Plots of 1000 simulations for each Scenario ########################
################################################################################################

## Make sure you choose the right LM for either MTN or WLG at the start of this section before you plot the graphs below:

#######################################################
############## PLOT FOR Scenarios A to I ##############
#######################################################

simObjects <- c("tempA", "tempB", "tempC", "tempD", "tempE", "tempF", "tempG", "tempH", "tempI")
detObjects <- c("N_projected_detA", "N_projected_detB", "N_projected_detC", "N_projected_detD", "N_projected_detE", "N_projected_detF", "N_projected_detG", "N_projected_detH", "N_projected_detI")

par(mfrow=c(3,3), oma=c(0,0,0,0), mar=c(5,4,2,1), las=1, bty="l")
maxY <- 120 ## max y-axis value
time <- 0:50 ## time interval for the plots

for(j in 1:9){
probExt <- prob_50years_wlg_lm[j,6]
#probExt <- prob_50years_mtn_lm[,6]##uncomment this line for mountain gorilla demographic parameter simulations
tempX <- get(simObjects[j])
N_projected_detX <- get(detObjects[j])
plot(N_projected_det~time, type="l", col=1, xlab="Years post-introduction", ylab="Population Size", ylim=c(0,maxY),lty=2, cex.lab=1, cex.axis=1, font.lab=2) ## plot of deterministic projection
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

########################################################################################
############################ PART 3: Individual-based model (IBM) ############################
########################################################################################

## For this model, we specified 5 categories of individuals: 
## immature (I): either weaned juveniles or unweaned infants, cycling adult female (C), pregnant adult female (P), lactating adult female (L), and females that lose their infants (CD)
## Thus, our IBM will inclde both an age and reproductive status category vector. See STEPS 1-2 below.
## Make sure to choose between either A. MTN or B. WLG for the steps below. 

########################################################################################
######################## STEP 1: AGE CATEGORIES FOR INDIVIDUALS ########################
########################################################################################

ages0 <- numeric(0) ## we need a vector that keeps track of the age of individuals
ages0 <- rep(1:length(n_mtn), rmultinom(1,100,prob=n_mtn)) ## We start with 100 individuals. rmultinom creates a vector that assigns the number of individuals (in this case 100) to each age, according to the LM demography (n or n_mtn)
hist(ages0)

## A. Parameters for optimistic model (MTN)
## Set up a variable that keeps track of the time since entry into the current age category
time0 <- numeric(length(ages0))

## Given an individual's age0 at time0, how many years has it been since entering the current age category? 
time0[(ages0-3.5)>=0 & (ages0-3.5)<4.5] <- ages0[(ages0-3.5)>=0 & (ages0-3.5)<4.5]-3.5 
time0[(ages0-8)>=0] <- 0
time0 <- time0[ages0>=3.5]
ages0 <- ages0[ages0>=3.5]

## B. Parameters for conservative model (WLG)
## Set up a variable that keeps track of the time since entry into the current age category
time0 <- numeric(length(ages0))

## Given an individual's age0 at time0, how many years has it been since entering the current age category? 
time0[(ages0-4.5)>=0 & (ages0-4.5)<5.5] <- ages0[(ages0-4.5)>=0 & (ages0-4.5)<5.5]-4.5 
time0[(ages0-10)>=0] <- 0
time0 <- time0[ages0>=4.5]
ages0 <- ages0[ages0>=4.5]

########################################################################################
######################## STEP 2: BREEDING STATUS OF INDIVIDUALS ########################
########################################################################################

## A. Adult female status (MTN)
status0 <- character(length(ages0))
status0[ages0<8] <- "I" ## Immature individuals are those between 0-8 years old and non-reproductive
status0[ages0>=8] <- "C" ## We start the simulation with cycling females only. With each timestep in the model, C females may transition to other categories, detailed in Step 3
data.frame(ages0, time0, status0)

## B. Adult female status (WLG)
status0 <- character(length(ages0))
status0[ages0<10] <- "I" ## Immature individuals are those between 0-10 years old and non-reproductive
status0[ages0>=10] <- "C" ## We start the simulation with cycling females only. With each timestep in the model, C females may transition to other categories, detailed in Step 3
data.frame(ages0, time0, status0)

## So, there are 3 vectors that we need to update at each time step: ages0, time0, status0 

########################################################################################
#################### STEP 3: TRANSITION PROBABILITES OF INDIVIDUALS ####################
########################################################################################

## Next, we specifiy the transition probabilities, in relation to time since entry in the "current" category

## A. Initial parameters: survivorship, fertility, and weaning age (MTN)
datX <- dat1 ## from the LM loaded in Part 2
weaningAge <- 4.5

## B. Initial parameters: survivorship, fertility, and weaning age (WLG)
#datX <- dat1 ## from the LM loaded in Part 2
#weaningAge <- 3.5

## Transition probabilities (MTN or WLG)
## t is the time passed in the initial class
timeunit <- 1/12 ## in years. Each time step is 1 month.

## 1. Transition from I to C:
IC <- function(t) ifelse(t<5.5, 0, 1) ## Probability for an subadult to transition to adult is zero if time as subadult is t<5.5 and 1 if t>=5.5

## 2. Transition from C to P:
alpha <- .99 ## Probability to be pregnant after 12 months being cycling
## relationship between alpha and p: 
## (1-p)^(1/timeunit)=1-alpha
## log(1-p)=log(1-alpha)/(1/timeunit)
## p=1-exp(log(1-alpha)/(1/timeunit))
CP <- function(t, alpha) 1-exp(log(1-alpha)/(1/timeunit)) ## Probability for a cycling adult female to become pregnant, per time step, given the alpha value

## 3. Transition from P to L:
PL <- function(t) ifelse(t<8.5/12, 0, 1) ## Probability for a pregnant adult female to transition to "lactating" is 1 if she has been pregnant for 8.5 months

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
alpha <- 0.85

simTshia <- function(ages0, status0, time0, nyears=50, timeunit=1/12, alpha=alpha, verbose=T){
  abmDataLog <- data.frame(timestep=0, ages=ages0, time=time0, status=status0, indiv=1:length(ages0), stringsAsFactors = FALSE)
  iter <- max(abmDataLog$indiv)+1 ## iter keeps track of the index of the last row (of abmDataLog) + 1, which is essentially the next available gorilla "ID" in the model. We don't want to repeat ID's. 
  for(i in 1:trunc(nyears/timeunit)){
    if(i%%(1/timeunit)==0 & verbose) print(paste("time =",i*timeunit)) ## %% or modulo find the remainder: every year, when you divide 1 by 12 %% = 0, so this keeps track of time in the simulation. if we select verbose=F, the time is not printed. 
    newAbmData <- abmDataLog[0,]
    abmData <- abmDataLog[abmDataLog$timestep==(i-1),] ## abmData keeps tracks of data from previous timestep
    for(j in 1:nrow(abmData)){ ## j is the individual in the model
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
      if(currentStatus=="L" & newStatus=="C" & infantSex==1) { ## 50/50 sex ratio
        newAbmData <- rbind(newAbmData, data.frame(timestep=i, ages=weaningAge+timeunit, time=weaningAge+timeunit, status="I", indiv=iter, stringsAsFactors = FALSE))
        iter <- iter+1
      } ## If weaning, add new row to abmData
    }
    if(sum(newAbmData$status=="CD")>0) newAbmData$status[newAbmData$status=="CD"] <- "C"
    if(nrow(newAbmData[newAbmData$status!="D",])==0) break ## if the adult female died, newAbmData should have 1 row with status D, and the juvenile gets weaned, so we keep the baby in the next time step. 
    abmDataLog <- rbind(abmDataLog, newAbmData[newAbmData$status!="D",]) ## if there are individuals that are not dead then add them to the newAbmData
  }
  return(abmDataLog)
}

########################################################################################
######################### STEP 5: RUN REINTRODUCTION SCENARIOS #########################
########################################################################################

nyears <- 20
timeunit <- 1/12
nruns <-1
alpha <- 0.85 ## alpha = 0.42 for WLG
res <- data.frame(Year=numeric(0), nyears=numeric(0), nInf = numeric(0))

for(i in 1:nruns){
  abmDataLog <- simTshia(ages0 = rep(10, 100), status0 =rep(c("C", "L"), c(20, 80)) , time0 = 0, nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=T)
  abmDataLog2 <- abmDataLog[abmDataLog$ages>=10,]
  for (j in unique(abmDataLog2$indiv)){
  	temp <- abmDataLog2[abmDataLog2$indiv==j & abmDataLog2$timestep>(10/timeunit),"status"]
  	nInf <- sum(((temp=="L")[-1]-(temp=="L")[-length(temp)])==1)
  	nYearObs <- length(temp)*timeunit
  	res <- rbind(res, data.frame(nyears= nYearObs, nInf = nInf))
  } 
}

apply(res, 2, sum)[2]/apply(res, 2, sum)[1]

## Now, let's apply the functions using each reintrodcution scenario. Load the appropriate csv file.
ReintroScenario_IBM <- read.csv(paste0(workingDir, "ReintroductionScenarios_IBM.csv"))

ReintroScenario_list <- list(ReA = list(Ages=ReintroScenario_IBM$ReA, Status=NA), ReB = list(Ages=ReintroScenario_IBM$ReB,  Status=NA), ReC = list(Ages=ReintroScenario_IBM$ReC, Status=NA), ReD = list(Ages=ReintroScenario_IBM$ReD, Status=NA), ReE = list(Ages=ReintroScenario_IBM$ReE, Status=NA), ReF = list(Ages=ReintroScenario_IBM$ReF, Status=NA), ReG = list(Ages=ReintroScenario_IBM$ReG, Status=NA), ReH = list(Ages=ReintroScenario_IBM$ReH, Status=NA), ReI = list(Ages=ReintroScenario_IBM$ReI, Status=NA)) ## create a list of these scenarios since columns are of different lengths
for(i in 1:length(ReintroScenario_list)){
  for(j in 1:length(ReintroScenario_list[[i]][j])){
    if(ReintroScenario_list[[1]][[2]][1] < 8) {
      ReintroScenario_list[[1]][[2]][1] <- "I" 
    } else {
      ReintroScenario_list[[1]][[2]][1] <- "C"
    }
  }
}


ReintroScenario_list <- lapply(ReintroScenario_list, function(x) x[!is.na(x)]) ## remove NAs

## Scenarios A1-I1: Mountain gorilla
## Scenarios A2-I2: Western gorilla

## Scenario A
## Scenario A1: one individual (adult female 19 years old, no reintroduced females added to the population)
## ages0 = c(19), status0 = c("C"), time0 = c(11), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
## Scenario A2: one individual (adult female 19 years old, no reintroduced females added to the population)
## ages0 = c(19), status0 = c("C"), time0 = c(9), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
nyears <- 50
timeunit <- 1/12
nruns <-100
alpha <- 0.99
res <- matrix(0, nrow=trunc(nyears/timeunit)+1, ncol=nruns)
for(i in 1:nruns){
  print(i)
  abmDataLog <- simTshia(ages0 = 19, status0 ="C", time0 = 0, nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F)
  nindiv <- tapply(abmDataLog$status,abmDataLog$timestep, function(v) length(v)+rbinom(1, sum(v=="L"), .5))##we're adding the unweaned females
  res[1:length(nindiv),i] <- nindiv
}
dim(res)

write.csv(res, file=paste0(workingDir,"pva_ABM_50year/ScenarioA.csv"), row.names=F)

## Scenario B
## Scenario B1: two individuals (adult female 19 years old, juvenile is a female, no reintroduced females added to the population)
## ages0 = c(5,19), status0 =c("I","C"), time0 = c(1.5, 11), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
## Scenario B2: two individual (adult female 19 years old, no reintroduced females added to the population)
## ages0 = c(5,19), status0 = c("I","C"), time0 = c(0.5, 9), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
nyears <- 50
timeunit <- 1/12
nruns <-1000
alpha=0.99
res0 <- matrix(0, nrow=trunc(nyears/timeunit)+1, ncol=nruns)
for(i in 1:nruns){
  print(i)
  abmDataLog <- simTshia(ages0 = c(5,19), status0 =c("I","C"), time0 = c(.5, 0), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F)
  nindiv <- tapply(abmDataLog$status,abmDataLog$timestep, function(v) length(v)+rbinom(1, sum(v=="L"), .5))##we're adding the unweaned females
  res0[1:length(nindiv),i] <- nindiv
}
res0

write.csv(res0, file=paste0(workingDir,"pva_ABM_50year/ScenarioB.csv"), row.names=F)

## Scenario C
## Scenario C1: three individuals (adult female 19 years old, 2 reintroduced females added to the population)
## ages0 = c(7,7,19), status0 = c("I","I","C"), time0 = c(3.5, 3.5, 11), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
## Scenario C2: three individuals (adult female 19 years old, 2 reintroduced females added to the population)
## ages0 = c(7,7,19), status0 = c("I","I","C"), time0 = c(2.5, 2.5, 9), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
nyears <- 50
timeunit <- 1/12
nruns <-1000
alpha=0.99
res1 <- matrix(0, nrow=trunc(nyears/timeunit)+1, ncol=nruns)
for(i in 1:nruns){
  print(i)
  abmDataLog <- simTshia(ages0 = c(7,7,19), status0 = c("I","I","C"), time0 = c(3.5, 3.5, 11), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F)
  nindiv <- tapply(abmDataLog$status,abmDataLog$timestep, function(v) length(v)+rbinom(1, sum(v=="L"), .5))##we're adding the unweaned females
  res1[1:length(nindiv),i] <- nindiv
}
res1

write.csv(res1, file=paste0(workingDir,"pva_ABM_50year/ScenarioC.csv"), row.names=F)

## Scenario D
## Scenario D1: four individuals (adult female 19 years old, 3 reintroduced females added to the population)
## ages0 = c(7,7,8,19), status0 = c("I","I","I","C"), time0 = c(3.5, 3.5, 4.5, 11), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
## Scenario D2: four individuals (adult female 19 years old, 3 reintroduced females added to the population)
## ages0 = c(7,7,8,19), status0 = c("I","I","I","C"), time0 = c(2.5, 2.5, 3.5, 9), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
nyears <- 50
timeunit <- 1/12
nruns <-1000
alpha=0.99
res2 <- matrix(0, nrow=trunc(nyears/timeunit)+1, ncol=nruns)
for(i in 1:nruns){
  print(i)
  abmDataLog <- simTshia(ages0 = c(7,7,8,19), status0 = c("I","I","I","C"), time0 = c(3.5, 3.5, 4.5, 11), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F)
  nindiv <- tapply(abmDataLog$status,abmDataLog$timestep, function(v) length(v)+rbinom(1, sum(v=="L"), .5))##we're adding the unweaned females
  res2[1:length(nindiv),i] <- nindiv
}
res2

write.csv(res2, file=paste0(workingDir,"pva_ABM_50year/ScenarioD.csv"), row.names=F)

## Scenario E
## Scenario E1: five individuals (adult female 19 years old, 4 reintroduced females added to the population)
## ages0 = c(7,7,8,9,19), status0 =c("I","I","I","C","C"), time0 = c(3.5, 3.5, 4.5, 5.5, 11), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
## Scenario E2: five individuals (adult female 19 years old, 4 reintroduced females added to the population)
## ages0 = c(7,7,8,9,19), status0 =c("I","I","I","C","C"), time0 = c(2.5, 2.5, 3.5, 4.5, 9), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
nyears <- 50
timeunit <- 1/12
nruns <-1000
alpha=0.99
res3 <- matrix(0, nrow=trunc(nyears/timeunit)+1, ncol=nruns)
for(i in 1:nruns){
  print(i)
  abmDataLog <- simTshia(ages0 = c(7,7,8,9,19), status0 =c("I","I","I","C","C"), time0 = c(3.5, 3.5, 4.5, 5.5, 11), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F)
  nindiv <- tapply(abmDataLog$status,abmDataLog$timestep, function(v) length(v)+rbinom(1, sum(v=="L"), .5))##we're adding the unweaned females
  res3[1:length(nindiv),i] <- nindiv
}
res3

write.csv(res3, file=paste0(workingDir,"pva_ABM_50year/ScenarioE.csv"), row.names=F)

## Scenario F
## Scenario F1: six individuals (adult female 19 years old, 5 reintroduced females added to the population)
## ages0 = c(7,7,8,9,12,19), status0 =c("I","I","I","C","C","C"), time0 = c(3.5, 3.5, 4.5, 5.5, 4, 11), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
## Scenario F2: six individuals (adult female 19 years old, 5 reintroduced females added to the population)
## ages0 = c(7,7,8,9,12,19), status0 =c("I","I","I","C","C","C"), time0 = c(2.5, 2.5, 3.5, 4.5, 2, 9), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
nyears <- 50
timeunit <- 1/12
nruns <-1000
alpha=0.99
res4 <- matrix(0, nrow=trunc(nyears/timeunit)+1, ncol=nruns)
for(i in 1:nruns){
  print(i)
  abmDataLog <- simTshia(ages0 = c(7,7,8,9,12,19), status0 =c("I","I","I","C","C","C"), time0 = c(3.5, 3.5, 4.5, 5.5, 4, 11), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F)
  nindiv <- tapply(abmDataLog$status,abmDataLog$timestep, function(v) length(v)+rbinom(1, sum(v=="L"), .5))##we're adding the unweaned females
  res4[1:length(nindiv),i] <- nindiv
}
res4

write.csv(res4, file=paste0(workingDir,"pva_ABM_50year/ScenarioF.csv"), row.names=F)

## Scenario G
## Scenario G1: seven individuals (adult female 19 years old, 6 reintroduced females added to the population)
## ages0 = c(7,7,8,9,12,12,19), status0 =c("I","I","I","C","C","C","C"), time0 = c(3.5, 3.5, 4.5, 5.5, 4, 4, 11), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
## Scenario G2: seven individuals (adult female 19 years old, 6 reintroduced females added to the population)
## ages0 = c(7,7,8,9,12,12,19), status0 =c("I","I","I","C","C","C","C"), time0 = c(2.5, 2.5, 3.5, 4.5, 2, 2, 9), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
nyears <- 50
timeunit <- 1/12
nruns <-1000
alpha=0.99
res5 <- matrix(0, nrow=trunc(nyears/timeunit)+1, ncol=nruns)
for(i in 1:nruns){
  print(i)
  abmDataLog <- simTshia(ages0 = c(7,7,8,9,12,12,19), status0 =c("I","I","I","C","C","C","C"), time0 = c(3.5, 3.5, 4.5, 5.5, 4, 4, 11), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F)
  nindiv <- tapply(abmDataLog$status,abmDataLog$timestep, function(v) length(v)+rbinom(1, sum(v=="L"), .5))##we're adding the unweaned females
  res5[1:length(nindiv),i] <- nindiv
}
res5

write.csv(res5, file=paste0(workingDir,"pva_ABM_50year/ScenarioG.csv"), row.names=F)

## Scenario H
## Scenario H1: eight individuals (adult female 19 years old, 7 reintroduced females added to the population)
## ages0 = c(7,7,8,9,12,12,17,19), status0 =c("I","I","I","C","C","C","C","C"), time0 = c(3.5, 3.5, 4.5, 5.5, 4, 4, 9, 11), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
## Scenario H2: eight individuals (adult female 19 years old, 7 reintroduced females added to the population)
## ages0 = c(7,7,8,9,12,12,17,19), status0 =c("I","I","I","C","C","C","C","C"), time0 = c(3.5, 3.5, 4.5, 5.5, 2, 2, 7, 9), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
nyears <- 50
timeunit <- 1/12
nruns <-1000
alpha=0.99
res6 <- matrix(0, nrow=trunc(nyears/timeunit)+1, ncol=nruns)
for(i in 1:nruns){
  print(i)
  abmDataLog <- simTshia(ages0 = c(7,7,8,9,12,12,17,19), status0 =c("I","I","I","C","C","C","C","C"), time0 = c(3.5, 3.5, 4.5, 5.5, 4, 4, 9, 11), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F)
  nindiv <- tapply(abmDataLog$status,abmDataLog$timestep, function(v) length(v)+rbinom(1, sum(v=="L"), .5))##we're adding the unweaned females
  res6[1:length(nindiv),i] <- nindiv
}
res6

write.csv(res6, file=paste0(workingDir,"pva_ABM_50year/ScenarioH.csv"), row.names=F)

## Scenario I
## Scenario I1: nine individuals (adult female 19 years old, 8 reintroduced females added to the population)
## ages0 = c(7,7,8,9,12,12,17,17,19), status0 =c("I","I","I","C","C","C","C","C","C"), time0 = c(3.5, 3.5, 4.5, 5.5, 4, 4, 9, 9, 11), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
## Scenario I2: nine individuals (adult female 19 years old, 8 reintroduced females added to the population)
## ages0 = c(7,7,8,9,12,12,17,17,19), status0 =c("I","I","I","C","C","C","C","C"), time0 = c(3.5, 3.5, 4.5, 5.5, 2, 2, 7, 7, 9), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F
nyears <- 50
timeunit <- 1/12
nruns <-1000
alpha=0.99
res7 <- matrix(0, nrow=trunc(nyears/timeunit)+1, ncol=nruns)
for(i in 1:nruns){
  print(i)
  abmDataLog <- simTshia(ages0 = c(7,7,8,9,12,12,17,17,19), status0 =c("I","I","I","C","C","C","C","C","C"), time0 = c(3.5, 3.5, 4.5, 5.5, 4, 4, 9, 9, 11), nyears=nyears, alpha=alpha, timeunit=timeunit, verbose=F)
  nindiv <- tapply(abmDataLog$status,abmDataLog$timestep, function(v) length(v)+rbinom(1, sum(v=="L"), .5))##we're adding the unweaned females
  res7[1:length(nindiv),i] <- nindiv
}
res7

write.csv(res7, file=paste0(workingDir,"pva_ABM_50year/ScenarioI.csv"), row.names=F)

########################################################################################
################# STEP 6: EXAMINE DATA FROM REINTRODUCTION SCENARIOS ###################
########################################################################################

### Now that the csv files have been written, we may not want to re-run the code for as long in the future, so we can just read the generated files and plot the data directly. Make sure you're reading the csv files from the correct folder. 

## Select the correct folder for either WLG or MTN data
workingDir_IBM <- "~/Box Sync/PVA_Paper/PVA_Tshiaberimu_R/pva_ABM_50year_Breuer_wlg/"
##workingDir_IBM <- "~/Box Sync/PVA_Paper/PVA_Tshiaberimu_R/pva_ABM_50year_Bronikowski_mtn/"

res <- as.matrix(read.csv(paste0(workingDir_IBM, "noIntro_1F.csv")))
finalPopSizes <- res[nrow(res),]

res0 <- as.matrix(read.csv(paste0(workingDir_IBM, "noIntro_2F.csv")))
finalPopSizes0 <- res0[nrow(res0),]

res1 <- as.matrix(read.csv(paste0(workingDir_IBM, "2Intro_3F.csv")))
finalPopSizes1 <- res1[nrow(res1),]

res2 <- as.matrix(read.csv(paste0(workingDir_IBM, "3Intro_4F.csv")))
finalPopSizes2 <- res2[nrow(res2),]

res3 <- as.matrix(read.csv(paste0(workingDir_IBM, "4Intro_5F.csv")))
finalPopSizes3 <- res3[nrow(res3),]

res4 <- as.matrix(read.csv(paste0(workingDir_IBM, "5Intro_6F.csv")))
finalPopSizes4 <- res4[nrow(res4),]

res5 <- as.matrix(read.csv(paste0(workingDir_IBM, "6Intro_7F.csv")))
finalPopSizes5 <- res5[nrow(res5),]

res6 <- as.matrix(read.csv(paste0(workingDir_IBM, "7Intro_8F.csv")))
finalPopSizes6 <- res6[nrow(res6),]

res7 <- as.matrix(read.csv(paste0(workingDir_IBM, "8Intro_9F.csv")))
finalPopSizes7 <- res7[nrow(res7),]

## what is probability that simulation results in extinction?
## see how many times final population is less than 0 (in last row of temp matrix or in all years of the projection?)

extA <- res[nrow(res),]==0
probExtA <- mean(res[nrow(res),]==0)
## OR probExtA <- length(extA[extA==TRUE])/length(extA)

extB <- res0[nrow(res0),]==0
probExtB <- mean(res0[nrow(res0),]==0)
## OR probExtB <- length(extB[extB==TRUE])/length(extB)

extC <- res1[nrow(res1),]==0
probExtC <- mean(res1[nrow(res1),]==0)
## OR probExtC <- length(extC[extC==TRUE])/length(extC)

extD <- res2[nrow(res2),]==0
probExtD <- mean(res2[nrow(res2),]==0)
## OR probExtD <- length(extD[extD==TRUE])/length(extD)

extE <- res3[nrow(res3),]==0
probExtE <- mean(res3[nrow(res3),]==0)
## OR probExtE <- length(extE[extE==TRUE])/length(extE)

extF <- res4[nrow(res4),]==0
probExtF <- mean(res4[nrow(res4),]==0)
## probExtF <- length(extF[extF==TRUE])/length(extF)

extG <- res5[nrow(res5),]==0
probExtG <- mean(res5[nrow(res5),]==0)
## probExtG <- length(extG[extG==TRUE])/length(extG)

extH <- res6[nrow(res6),]==0
probExtH <- mean(res6[nrow(res6),]==0)
## OR probExtH <- length(extH[extH==TRUE])/length(extH)

extI <- res7[nrow(res7),]==0
probExtI <- mean(res7[nrow(res7),]==0)
## OR probExtI <- length(extI[extI==TRUE])/length(extI)

##############################################################
############## EXAMINE LAMBDA VALUES FOR GROWTH ##############
##############################################################

## Scenario A
logLambda <- (1/50)*log(res[nrow(res),]/res[1,]) ## 50 years for the projected time period, lambda = 1/timeperiod*log(Ntfinal)/Nt0
lambda <- exp(logLambda)
lambdaA <- mean(lambda)

## Scenario B
logLambda0 <- (1/50)*log(res0[nrow(res0),]/res0[1,]) ## 50 years for the projected time period, lambda = 1/timeperiod*log(Ntfinal)/Nt0
lambda0 <- exp(logLambda0)
lambdaB <- mean(lambda0)

## Scenario C
logLambda1 <- (1/50)*log(res1[nrow(res1),]/res1[1,]) ## 50 years for the projected time period, lambda = 1/timeperiod*log(Ntfinal)/Nt0
lambda1 <- exp(logLambda1)
lambdaC <- mean(lambda1)

## Scenario D
logLambda2 <- (1/50)*log(res2[nrow(res2),]/res2[1,]) ## 50 years for the projected time period, lambda = 1/timeperiod*log(Ntfinal)/Nt0
lambda2 <- exp(logLambda2)
lambdaD <- mean(lambda2)

## Scenario E
logLambda3 <- (1/50)*log(res3[nrow(res3),]/res3[1,]) ## 50 years for the projected time period, lambda = 1/timeperiod*log(Ntfinal)/Nt0
lambda3 <- exp(logLambda3)
lambdaE <- mean(lambda3)

## Scenario F
logLambda4 <- (1/50)*log(res4[nrow(res4),]/res4[1,]) ## 50 years for the projected time period, lambda = 1/timeperiod*log(Ntfinal)/Nt0
lambda4 <- exp(logLambda4)
lambdaF <- mean(lambda4)

## Scenario G
logLambda5 <- (1/50)*log(res5[nrow(res5),]/res5[1,]) ## 50 years for the projected time period, lambda = 1/timeperiod*log(Ntfinal)/Nt0
lambda5 <- exp(logLambda5)
lambdaG <- mean(lambda5)

## Scenario H
logLambda6 <- (1/50)*log(res6[nrow(res6),]/res6[1,]) ## 50 years for the projected time period, lambda = 1/timeperiod*log(Ntfinal)/Nt0
lambda6 <- exp(logLambda6)
lambdaH <- mean(lambda6)

## Scenario I
logLambda7 <- (1/50)*log(res7[nrow(res7),]/res7[1,]) ## 50 years for the projected time period, lambda = 1/timeperiod*log(Ntfinal)/Nt0
lambda7 <- exp(logLambda7)
lambdaI <- mean(lambda7)

## let's save these lambda and prob extn values as a data frame and save it as a csv file
lambda_extn_ibm <- data.frame(scenario = as.factor(LETTERS[1:9]), lambda_50years =  c(lambdaA, lambdaB, lambdaC, lambdaD, lambdaE, lambdaF, lambdaG, lambdaH, lambdaI), prob_Extn =  (c(probExtA, probExtB, probExtC, probExtD, probExtE, probExtF, probExtG, probExtH, probExtI)*100))

## Make sure you use the correct LM for WLG or MTN when you run the code
## WLG values
write.csv(lambda_extn_ibm, file="/Users/neethaiyer/Box Sync/PVA_Paper/PVA_Tshiaberimu_R/pva_lambda_extn/lambda_extn_ibm_WLG.csv", row.names=F)
## MTN values
write.csv(lambda_extn_ibm, file="/Users/neethaiyer/Box Sync/PVA_Paper/PVA_Tshiaberimu_R/pva_lambda_extn/lambda_extn_ibm_MTN.csv", row.names=F)

################################################################################################
###################################### PLOTS for IBM ###########################################
################################################################################################
## You can save the plots using the commented out lines of code below, or directly by saving the plots using RStudio's interface. 
## file_path <- file.path("~.png",paste("PVA",i, ".png", sep = ""))
## png(file_path, width=4, height=4, units='in', res=100)

#################################################
########## PLOT FOR mean projections ############
#################################################

timeunit <- 1/12

plot((0:(nrow(res)-1))*timeunit, apply(res, 1, mean), type="l", col="azure4", lwd=3, ylim=c(0,120))
lines((0:(nrow(res0)-1))*timeunit, apply(res0, 1, mean), type="l", col="azure4", lwd=3)
lines((0:(nrow(res1)-1))*timeunit, apply(res1, 1, mean), type="l", col="azure4", lwd=3)
lines((0:(nrow(res2)-1))*timeunit, apply(res2, 1, mean), type="l", col="azure4", lwd=3)
lines((0:(nrow(res3)-1))*timeunit, apply(res3, 1, mean), type="l", col="azure4", lwd=3)
lines((0:(nrow(res4)-1))*timeunit, apply(res4, 1, mean), type="l", col="azure4", lwd=3)
lines((0:(nrow(res5)-1))*timeunit, apply(res5, 1, mean), type="l", col="azure4", lwd=3)
lines((0:(nrow(res6)-1))*timeunit, apply(res6, 1, mean), type="l", col="azure4", lwd=3)
lines((0:(nrow(res7)-1))*timeunit, apply(res7, 1, mean), type="l", col="azure4", lwd=3)

#################################################
######### PLOT of N final per Scenario ##########
#################################################
## extract the N finals of the mean population projections
s <- tail(apply(res, 1, mean), n=1)
s0 <- tail(apply(res0, 1, mean), n=1)
s1 <- tail(apply(res1, 1, mean), n=1)
s2 <- tail(apply(res2, 1, mean), n=1)
s3 <- tail(apply(res3, 1, mean), n=1)
s4 <- tail(apply(res4, 1, mean), n=1)
s5 <- tail(apply(res5, 1, mean), n=1)
s6 <- tail(apply(res6, 1, mean), n=1)
s7 <- tail(apply(res7, 1, mean), n=1)

## Extract upper and lower 95 intervals for population projections
c <- tail(apply(res, 1, function(v) quantile(v, probs=c(0.05, 0.95))), n=1)
test <- apply(res, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
apply(res0, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
apply(res1, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
apply(res2, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
apply(res3, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
apply(res4, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
apply(res5, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
apply(res6, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
apply(res7, 1, function(v) quantile(v, probs=c(0.05, 0.95)))

## Let's save those values as a dataframe
Nfinal_ibm <- data.frame(scenario = as.factor(LETTERS[1:9]), Nfinal_mean = c(s, s0, s1, s2, s3, s4, s5, s6, s7))

##### Plot of final population sizes for WLG
plot.window(xlim=c(1,9), ylim=c(0,150))
boxplot(res[nrow(res),], res0[nrow(res0),], res1[nrow(res1),], res2[nrow(res2),], res3[nrow(res3),], res4[nrow(res4),], res5[nrow(res5),], res6[nrow(res6),], res7[nrow(res7),], pch=21, col="coral", bg=coral, ylim=c(0,150), varwidth=FALSE)
lines(x=c(-5:51), y=rep(50, 57), col="red", lwd=2, lty=2)
par(new = TRUE, mar=c(5.1,4.1,4.1,5.1))
axis(2, at=seq(0, 150, by=50), labels=seq(0, 150, by=50), ylab="Population size after 50 years", font.lab=2)
axis(1, LETTERS[1:9], at=1:9, labels=LETTERS[1:9])
title(xlab="Reintroduction Scenario", ylab="Population size after 50 years", font.lab=2)

##### Plot of final population sizes for MTN
plot.window(xlim=c(1,9), ylim=c(0,150))
boxplot(res[nrow(res),], res0[nrow(res0),], res1[nrow(res1),], res2[nrow(res2),], res3[nrow(res3),], res4[nrow(res4),], res5[nrow(res5),], res6[nrow(res6),], res7[nrow(res7),], pch=24, col="azure4", bg=azure4, ylim=c(0,150), varwidth=FALSE)
lines(x=c(-5:51), y=rep(50, 57), col="red", lwd=2, lty=2)
par(new = TRUE, mar=c(5.1,4.1,4.1,5.1))
axis(2, at=seq(0, 150, by=50), labels=seq(0, 150, by=50), ylab="Population size after 50 years", font.lab=2)
axis(1, LETTERS[1:9], at=1:9, labels=LETTERS[1:9])
title(xlab="Reintroduction Scenario", ylab="Population size after 50 years", font.lab=2)

## Let's plot lambda and extinction probabilities
parameters_wlg_ibm <- as.data.frame(read.csv("/Users/neethaiyer/Box Sync/PVA_Paper/PVA_Tshiaberimu_R/pva_lambda_extn/lambda_extn_ibm_WLG.csv"))
parameters_mtn_ibm <- as.data.frame(read.csv("/Users/neethaiyer/Box Sync/PVA_Paper/PVA_Tshiaberimu_R/pva_lambda_extn/lambda_extn_ibm_MTN.csv"))

################################################################################################
################################ Plots of extinction probability  ##############################
################################# For WLG, MTN using LM and IBM  ###############################
################################################################################################
plot.new()
par(mar=c(5.1,4.1,4.1,5.1))
plot.window(xlim=c(1,9), ylim=c(0,100))
axis(1, 1:9, LETTERS[1:9])
axis(2)
axis(2, font.lab=2, at=seq(0, 100, by=10), labels=seq(0, 100, by=10))
title(xlab="Reintroduction Scenario", ylab="Probability of Extinction", font.lab=2)
lines(parameters_wlg_ibm$scenario, parameters_wlg_ibm$prob_Extn, bg="coral", type="b", pch=21)
lines(parameters_mtn_ibm$scenario, parameters_mtn_ibm$prob_Extn, bg="azure4", type="b", pch=24)
lines(parameters_wlg_lm$scenario, parameters_wlg_lm$prob_Extn, bg=coral, type="b", pch=21, lty=2)
lines(parameters_mtn_lm$scenario, parameters_mtn_lm$prob_Extn, bg=azure4, type="b", pch=24, lty=2)
legend(1, 100, legend=c("Western Gorillas - IBM", "Western Gorillas - LM", "Mountain Gorillas - IBM", "Mountain Gorillas - LM"),
       pt.bg=c("coral", coral, "azure4", azure4), pch=c(21,21,24,24), lty=c(1,2,1,2), cex=0.8, text.font=2)

################################################################################################
#################################### Plots of lambda values  ###################################
################################################################################################
plot.new()
par(mar=c(5.1,4.1,4.1,5.1))
plot.window(xlim=c(1,9), ylim=c(0,1.5))
axis(1, 1:9, LETTERS[1:9])
axis(2)
##axis(2, font.lab=2, at=seq(0, 1, by=0.1), labels=seq(0, 1, by=0.1))
title(xlab="Reintroduction Scenario", ylab="Lambda", font.lab=2)
lines(parameters_wlg_ibm$scenario, parameters_wlg_ibm$lambda_50years, bg="coral", type="b", pch=21)
lines(parameters_mtn_ibm$scenario, parameters_mtn_ibm$lambda_50years, bg="azure4", type="b", pch=24)
##lines(parameters_wlg_lm$scenario, parameters_wlg_lm$lambda_50years, bg="coral", type="b", pch=20, lty=2)
##lines(parameters_mtn_lm$scenario, parameters_mtn_lm$lambda_50years, bg="azure4", type="b", pch=20, lty=2)
legend(1, 1.4, legend=c("Western Gorillas - IBM", "Mountain Gorillas - IBM"),
       pt.bg=c("coral", "azure4"), pch=c(21,24), lty=c(1,1), cex=0.8, text.font=2)

################################################################################################
########################### Plots of 1000 simulations for each Scenario ########################
################################################################################################

## Make sure you choose the right workingDir_IBM for either WLG or MTN at the start of PART 6 section !!!!!
timeunit <- 1/12

#################################################
############## PLOT FOR Scenario A ##############
#################################################
par(mfrow=c(3,3))
plot((0:(nrow(res)-1))*timeunit, apply(res, 1, mean), type="l", col=2, lwd=2, xlab="Years post-introduction", ylab="Population size", font.lab=2, bty="l", ylim=c(0,120))
for(i in 1:ncol(res)){
  lines((0:(nrow(res)-1))*timeunit, res[,i], col=grey(.9))
}

## add mean trend
lines((0:(nrow(res)-1))*timeunit, apply(res, 1, mean), type="l", col=2, lwd=3)

## add 95% upper/lower limits
qtiles <- apply(res, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
lines((0:(nrow(res)-1))*timeunit, qtiles[1,], lty=2)
lines((0:(nrow(res)-1))*timeunit, qtiles[2,], lty=2)

## add some graph fluff
title(main="A: No introduction, male juvenile", sub=paste0("Probability of extinction = ", round(probExtA, 3)*100, "%"), cex.main=1, cex.sub=1, col.sub=1,  font.sub=3)

## add a line for the 50 individual mark
lines(x=c(-5:50), y=rep(50, 56), col="navy", lwd=2, lty=1)

#################################################
############## PLOT FOR Scenario B ##############
#################################################
plot((0:(nrow(res0)-1))*timeunit, apply(res0, 1, mean), type="l", col=2, lwd=2, xlab="Years post-introduction", ylab="Population size", font.lab=2, bty="l", ylim=c(0,120))
for(i in 1:ncol(res0)){
  lines((0:(nrow(res0)-1))*timeunit, res0[,i], col=grey(.9))
}

## add mean trend
lines((0:(nrow(res0)-1))*timeunit, apply(res0, 1, mean), type="l", col=2, lwd=3)

## add 95% upper/lower limits
qtiles <- apply(res0, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
lines((0:(nrow(res0)-1))*timeunit, qtiles[1,], lty=2)
lines((0:(nrow(res0)-1))*timeunit, qtiles[2,], lty=2)

## add some graph fluff
title(main="B: No introduction, female juvenile", sub=paste0("Probability of extinction = ", round(probExtB, 3)*100, "%"), cex.main=1, cex.sub=1, col.sub=1,  font.sub=3)

## add a line for the 50 individual mark
lines(x=c(-5:50), y=rep(50, 56), col="navy", lwd=2, lty=1)

#################################################
############## PLOT FOR Scenario C ##############
#################################################
plot((0:(nrow(res1)-1))*timeunit, apply(res1, 1, mean), type="l", col=2, lwd=2, xlab="Years post-introduction", ylab="Population size", font.lab=2, bty="l", ylim=c(0,120))
for(i in 1:ncol(res1)){
  lines((0:(nrow(res1)-1))*timeunit, res1[,i], col=grey(.9))
}

## add mean trend
lines((0:(nrow(res1)-1))*timeunit, apply(res1, 1, mean), type="l", col=2, lwd=3)

## add 95% upper/lower limits
qtiles <- apply(res1, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
lines((0:(nrow(res1)-1))*timeunit, qtiles[1,], lty=2)
lines((0:(nrow(res1)-1))*timeunit, qtiles[2,], lty=2)

## add some graph fluff
title(main="C: 2 females introduced", sub=paste0("Probability of extinction = ", round(probExtC, 3)*100, "%"), cex.main=1, cex.sub=1, col.sub=1,  font.sub=3)

## add a line for the 50 individual mark
lines(x=c(-5:50), y=rep(50, 56), col="navy", lwd=2, lty=1)

#################################################
############## PLOT FOR Scenario D ##############
#################################################
plot((0:(nrow(res2)-1))*timeunit, apply(res2, 1, mean), type="l", col=2, lwd=2, xlab="Years post-introduction", ylab="Population size", bty="l", ylim=c(0,120))
for(i in 1:ncol(res2)){
  lines((0:(nrow(res2)-1))*timeunit, res2[,i], col=grey(.9))
}

## add mean trend
lines((0:(nrow(res2)-1))*timeunit, apply(res2, 1, mean), type="l", col=2, lwd=3)

## add 95% upper/lower limits
qtiles <- apply(res2, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
lines((0:(nrow(res2)-1))*timeunit, qtiles[1,], lty=2)
lines((0:(nrow(res2)-1))*timeunit, qtiles[2,], lty=2)

## add some graph fluff
title(main="D: 3 females introduced", sub=paste0("Probability of extinction = ", round(probExtD, 3)*100, "%"), cex.main=1, cex.sub=1, col.sub=1,  font.sub=3)

## add a line for the 50 individual mark
lines(x=c(-5:50), y=rep(50, 56), col="navy", lwd=2, lty=1)

#################################################
############## PLOT FOR Scenario E ##############
#################################################
plot((0:(nrow(res3)-1))*timeunit, apply(res3, 1, mean), type="l", col=2, lwd=2, xlab="Years post-introduction", ylab="Population size", bty="l", ylim=c(0,120))
for(i in 1:ncol(res3)){
  lines((0:(nrow(res3)-1))*timeunit, res3[,i], col=grey(.9))
}

## add mean trend
lines((0:(nrow(res3)-1))*timeunit, apply(res3, 1, mean), type="l", col=2, lwd=3)

## add 95% upper/lower limits
qtiles <- apply(res3, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
lines((0:(nrow(res3)-1))*timeunit, qtiles[1,], lty=2)
lines((0:(nrow(res3)-1))*timeunit, qtiles[2,], lty=2)

## add some graph fluff
title(main="E: 4 females introduced", sub=paste0("Probability of extinction = ", round(probExtE, 3)*100, "%"), cex.main=1, cex.sub=1, col.sub=1,  font.sub=3)

## add a line for the 50 individual mark
lines(x=c(-5:50), y=rep(50, 56), col="navy", lwd=2, lty=1)

#################################################
############## PLOT FOR Scenario F ##############
#################################################
plot((0:(nrow(res4)-1))*timeunit, apply(res4, 1, mean), type="l", col=2, lwd=2, xlab="Years post-introduction", ylab="Population size", bty="l", ylim=c(0,120))
for(i in 1:ncol(res4)){
  lines((0:(nrow(res4)-1))*timeunit, res4[,i], col=grey(.9))
}

## add mean trend
lines((0:(nrow(res4)-1))*timeunit, apply(res4, 1, mean), type="l", col=2, lwd=3)

## add 95% upper/lower limits
qtiles <- apply(res4, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
lines((0:(nrow(res4)-1))*timeunit, qtiles[1,], lty=2)
lines((0:(nrow(res4)-1))*timeunit, qtiles[2,], lty=2)

## add some graph fluff
title(main="F: 5 females introduced", sub=paste0("Probability of extinction = ", round(probExtF, 3)*100, "%"), cex.main=1, cex.sub=1, col.sub=1,  font.sub=3)

## add a line for the 50 individual mark
lines(x=c(-5:50), y=rep(50, 56), col="navy", lwd=2, lty=1)

#################################################
############## PLOT FOR Scenario G ##############
#################################################
plot((0:(nrow(res5)-1))*timeunit, apply(res5, 1, mean), type="l", col=2, lwd=2, xlab="Years post-introduction", ylab="Population size", bty="l", ylim=c(0,120))
for(i in 1:ncol(res5)){
  lines((0:(nrow(res5)-1))*timeunit, res5[,i], col=grey(.9))
}

## add mean trend
lines((0:(nrow(res5)-1))*timeunit, apply(res5, 1, mean), type="l", col=2, lwd=3)

## add 95% upper/lower limits
qtiles <- apply(res5, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
lines((0:(nrow(res5)-1))*timeunit, qtiles[1,], lty=2)
lines((0:(nrow(res5)-1))*timeunit, qtiles[2,], lty=2)

## add some graph fluff
title(main="G: 6 females introduced", sub=paste0("Probability of extinction = ", round(probExtG, 3)*100, "%"), cex.main=1, cex.sub=1, col.sub=1,  font.sub=3)

## add a line for the 50 individual mark
lines(x=c(-5:50), y=rep(50, 56), col="navy", lwd=2, lty=1)

#################################################
############## PLOT FOR Scenario H ##############
#################################################
plot((0:(nrow(res6)-1))*timeunit, apply(res6, 1, mean), type="l", col=2, lwd=2, xlab="Years post-introduction", ylab="Population size", bty="l", ylim=c(0,120))
for(i in 1:ncol(res6)){
  lines((0:(nrow(res6)-1))*timeunit, res6[,i], col=grey(.9))
}

## add mean trend
lines((0:(nrow(res6)-1))*timeunit, apply(res6, 1, mean), type="l", col=2, lwd=3)

## add 95% upper/lower limits
qtiles <- apply(res6, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
lines((0:(nrow(res6)-1))*timeunit, qtiles[1,], lty=2)
lines((0:(nrow(res6)-1))*timeunit, qtiles[2,], lty=2)

## add some graph fluff
title(main="H: 7 females introduced", sub=paste0("Probability of extinction = ", round(probExtH, 3)*100, "%"), cex.main=1, cex.sub=1, col.sub=1,  font.sub=3)

## add a line for the 50 individual mark
lines(x=c(-5:50), y=rep(50, 56), col="navy", lwd=2, lty=1)

#################################################
############## PLOT FOR Scenario I ##############
#################################################
plot((0:(nrow(res7)-1))*timeunit, apply(res7, 1, mean), type="l", col=2, lwd=2, xlab="Years post-introduction", ylab="Population size", bty="l", ylim=c(0,120))
for(i in 1:ncol(res7)){
  lines((0:(nrow(res7)-1))*timeunit, res7[,i], col=grey(.9))
}

## add mean trend
lines((0:(nrow(res7)-1))*timeunit, apply(res7, 1, mean), type="l", col=2, lwd=3)

## add 95% upper/lower limits
qtiles <- apply(res7, 1, function(v) quantile(v, probs=c(0.05, 0.95)))
lines((0:(nrow(res7)-1))*timeunit, qtiles[1,], lty=2)
lines((0:(nrow(res7)-1))*timeunit, qtiles[2,], lty=2)

## add some graph fluff
title(main="I: 8 females introduced", sub=paste0("Probability of extinction = ", round(probExtI, 3)*100, "%"), cex.main=1, cex.sub=1, col.sub=1,  font.sub=3)

## add a line for the 50 individual mark
lines(x=c(-5:50), y=rep(50, 56), col="navy", lwd=2, lty=1)

################################################################################################
################################# Additional Exploratory plots  ################################
################################################################################################

## look at distribution of final population sizes for 1000 simulations
par(mfrow=c(3,3))
plot(table(finalPopSizes))
plot(table(finalPopSizes0))
plot(table(finalPopSizes1))
plot(table(finalPopSizes2))
plot(table(finalPopSizes3))
plot(table(finalPopSizes4))
plot(table(finalPopSizes5))
plot(table(finalPopSizes6))
plot(table(finalPopSizes7))

## look at distribution of final population sizes for 1000 simulations
par(mfrow=c(3,3))
hist(finalPopSizes)
abline(v=mean(finalPopSizes), col="red")
hist(finalPopSizes0)
abline(v=mean(finalPopSizes0), col="red")
hist(finalPopSizes1)
abline(v=mean(finalPopSizes1), col="red")
hist(finalPopSizes2)
abline(v=mean(finalPopSizes2), col="red")
hist(finalPopSizes3)
abline(v=mean(finalPopSizes3), col="red")
hist(finalPopSizes4)
abline(v=mean(finalPopSizes4), col="red")
hist(finalPopSizes5)
abline(v=mean(finalPopSizes5), col="red")
hist(finalPopSizes6)
abline(v=mean(finalPopSizes6), col="red")
hist(finalPopSizes7)
abline(v=mean(finalPopSizes7), col="red")
