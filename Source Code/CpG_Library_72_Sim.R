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
newmetric72max <- 0.01
#probabilities of a cpg being slected which will be updated at every iteration
#we start with equal probability of each cpg being selected
probs72 <- rep(1/length(candobject$candidateDMRs), length(candobject$candidateDMRs))
names(probs72) <- candobject$candidateDMRs
#register cores to run parts of the simulation in parallel 
registerDoParallel(16)
for(i in 1:nsim){
  #get libraries of fixed size and update probabilities
  #72
  random_cand72 <- update_probabilities_par(cpg_list=candobject$candidateDMRs, p=probs72,
                                            referencedata=RefBetas, cellident=cell_ident, fixed=T, n=72, ncores=16)
  
  library72 <- random_cand72[[1]]
  p72 <- random_cand72[[2]]
  #update the probabilities at the correct indexes
  probs72[library72] <- p72
  #scale the probabilities to sum to 1
  probs72 <- probs72/sum(probs72)
  #store the new metric
  newmetric72 <- random_cand72[[3]]
  
  #check to see if the new metric is larger than our previous max for each size library
  if(newmetric72 > newmetric72max){
    save(newmetric72, library72, file = "Optimal set of 72 CpGs.RData")
    newmetric72max <- newmetric72
    print(paste("The max new metric for library size 72 so far is: ", round(newmetric72max, 3), sep = ""))
  }
  
  #print an indicator of every 1000 iterations
  if(i %% 1000 == 0){
    print(paste("Iteration number: ", i, sep = ""))
  }
  
}