#####################################################################################
# install and load necessary R pacakges
#####################################################################################
install.packages(c("quadprog", "doParallel", "foreach", "genefilter"))
library(quadprog)
library(genefilter)
library(doParallel)
library(foreach)

#####################################################################################
# Read in necessary objects and data 
#####################################################################################
load("Reference Data for IDOL.RData")  # reference data
source("UnsupervisedOptimalDMRFinderV3.R") #new metric functions 


uniquecelltypes <- c("CD4T", "CD8T", "NK", "Bcell", "Mono", "Gran")
keep <- which(referencePd$CellType %in% uniquecelltypes)
referencePd.sub <- referencePd[keep,]
RefBetas <- referenceBetas[,keep]


#cell identities for the columns of refbetas dsc
cols <- colnames(RefBetas)
cell_ident <- sapply(strsplit(cols, split="_", fixed = T), function(x) (x[1]))


#####################################################################################
# Simulation for finding the best library with new metric
#####################################################################################

#note this candidate object will be used to test our metric
candobject = CandidateDMRFinder.v2(uniquecelltypes, referenceBetas, referencePd, M = 150)
#number of simulations to run
nsim <- 100000
#initialize to store max new metric result for simulation
newmetricrandommax <- 0.01
#probabilities of a cpg being slected which will be updated at every iteration
#we start with equal probability of each cpg being selected
probsrandom <- rep(1/length(candobject$candidateDMRs), length(candobject$candidateDMRs))
names(probsrandom) <- candobject$candidateDMRs

#register cores to run parts of the simulation in parallel 
registerDoParallel(16)
for(i in 1:nsim){
  #get libraries of random size and update probabilities
  random_candrand <- update_probabilities_par(cpg_list=candobject$candidateDMRs, p=probsrandom,
                                              referencedata=RefBetas, cellident=cell_ident, fixed=F, 
                                              lower=50,upper=1500, ncores=16)
  
  libraryrand <- random_candrand[[1]]
  prand <- random_candrand[[2]]
  #update the probabilities for the cpgs in the random library
  #here we are updating using the cpg names as the indexes
  probsrandom[libraryrand] <- prand
  #scale the probabilities to sum to 1
  probsrandom <- probsrandom/sum(probsrandom)
  newmetricrandom <- random_candrand[[3]]
  
  if(newmetricrandom > newmetricrandommax){
    save(newmetricrandom, libraryrand, file = "Optimal set of Random Number of CpGs.RData")
    newmetricrandommax <- newmetricrandom
    print(paste("The size of the best library so far is: ", length(libraryrand), sep = ""))
    print(paste("The max new metric for random library size so far is: ", round(newmetricrandommax, 3), sep = ""))
  }
  
  
  #print an indicator of every 1000 iterations
  if(i %% 1000 == 0){
    print(paste("Iteration number: ", i, sep = ""))
  }
  
}