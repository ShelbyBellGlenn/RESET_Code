#####################################################################################
# install and load necessary R pacakges
#####################################################################################
install.packages(c("quadprog", "genefilter", "doParallel", "foreach"))
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
newmetric540max <- 0.01
#probabilities of a cpg being slected which will be updated at every iteration
#we start with equal probability of each cpg being selected
probs540 <- rep(1/length(candobject$candidateDMRs), length(candobject$candidateDMRs))
names(probs540) <- candobject$candidateDMRs
#register cores to run parts of the simulation in parallel 
registerDoParallel(16)
for(i in 1:nsim){
  #540
  random_cand540 <- update_probabilities_par(cpg_list=candobject$candidateDMRs, p=probs540,
                                             referencedata=RefBetas, cellident=cell_ident, fixed=T, n=540, ncores=16)
  library540 <- random_cand540[[1]]
  p540 <- random_cand540[[2]]
  #update the probabilities at the correct indexes
  probs540[library540] <- p540
  #scale the probabilities to sum to 1
  probs540 <- probs540/sum(probs540)
  newmetric540 <- random_cand540[[3]]
  
  if(newmetric540 > newmetric540max){
    save(newmetric540, library540, file = "Optimal set of 540 CpGs.RData")
    newmetric540max <- newmetric540
    print(paste("The max new metric for library size 540 so far is: ", round(newmetric540max, 3), sep = ""))
  }
  
  #print an indicator of every 1000 iterations
  if(i %% 1000 == 0){
    print(paste("Iteration number: ", i, sep = ""))
  }
}