## Set the working directory:
workingDir <- "~/Box Sync/PVA_Paper/PVA_Tshiaberimu_R/"
##workingDir <- "~/Documents/git repositories/PVA_Tshiaberimu_R/"

############################################################################################
################################# LESLIE MATRIX PROJECTIONS ################################
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