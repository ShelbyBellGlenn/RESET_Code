#########Simulation results##########

#Load in training data
load("Reference Data for IDOL.RData")  # reference data
load("AdultMixed.RData")
load("MethodA.RData")
load("MethodB.RData")
source("UnsupervisedOptimalDMRFinderV3.R") #new metric functions 
source("IDOL revised functions 02_19_2017.R") #idol functions

#combine the three datasets
combinedbetas <- cbind(AdultMixed.betas, MethodA.betas, MethodB.betas)

#combine the training covariates
truepropscombined <- rbind(AdultMixed.covariates, MethodA.covariates, MethodB.covariates)


#get mean methylation values for each cell type from the orgininal ref betas matrix
uniquecelltypes <- c("CD4T", "CD8T", "NK", "Bcell", "Mono", "Gran")
keep <- which(referencePd$CellType %in% uniquecelltypes)
referencePd.sub <- referencePd[keep,]
RefBetas <- referenceBetas[,keep]
meanmeth = matrix(NA, nrow = 485512, ncol = 6)
rownames(meanmeth) = rownames(RefBetas)
colnames(meanmeth) = uniquecelltypes
for(k in 1:6) {
  ind = which(referencePd.sub$CellType %in% uniquecelltypes[k])
  meanmeth[,k] = apply(RefBetas[,ind], 1, mean, na.rm = T)
}

#cell identities for the columns of refbetas dsc
cols <- colnames(RefBetas)
cell_ident <- sapply(strsplit(cols, split="_", fixed = T), function(x) (x[1]))

#create RMSE function
rmse <- function(p, o){
  sqrt(sum((p - o)^2/length(p)))
}

#load our optimal libraries
load("Optimal set of 72 CpGs.RData")
load("Optimal set of 120 CpGs.RData")
load("Optimal set of 180 CpGs.RData")
load("Optimal set of 240 CpGs.RData")
load("Optimal set of 300 CpGs.RData")
load("Optimal set of 360 CpGs.RData")
load("Optimal set of 540 CpGs.RData")
load("Optimal set of Random Number of CpGs.RData")

#use full training data set and subset into the randomly selected cpgs
#for the new metric
training_sub72 <- combinedbetas[library72,]
training_sub120 <- combinedbetas[library120,]
training_sub180 <- combinedbetas[library180,]
training_sub240 <- combinedbetas[library240,]
training_sub300 <- combinedbetas[library300,]
training_sub360 <- combinedbetas[library360,]
training_sub540 <- combinedbetas[library540,]
training_sub_rand <- combinedbetas[libraryrand,]

#get the mean meth values
#for the new metric
meanmeth_sub72 <- meanmeth[library72,]
meanmeth_sub120 <- meanmeth[library120,]
meanmeth_sub180 <- meanmeth[library180,]
meanmeth_sub240 <- meanmeth[library240,]
meanmeth_sub300 <- meanmeth[library300,]
meanmeth_sub360 <- meanmeth[library360,]
meanmeth_sub540 <- meanmeth[library540,]
meanmeth_subrand <- meanmeth[libraryrand,]

#get predicted cell proportions for each cell types using training data
#for the new metric
Cell_estimates72 <- projectWBCnew(Y= training_sub72, coefWBC = meanmeth_sub72)*100
Cell_estimates120 <- projectWBCnew(Y= training_sub120, coefWBC = meanmeth_sub120)*100
Cell_estimates180 <- projectWBCnew(Y= training_sub180, coefWBC = meanmeth_sub180)*100
Cell_estimates240 <- projectWBCnew(Y= training_sub240, coefWBC = meanmeth_sub240)*100
Cell_estimates300 <- projectWBCnew(Y= training_sub300, coefWBC = meanmeth_sub300)*100
Cell_estimates360 <- projectWBCnew(Y= training_sub360, coefWBC = meanmeth_sub360)*100
Cell_estimates540 <- projectWBCnew(Y= training_sub540, coefWBC = meanmeth_sub540)*100
Cell_estimatesrand <- projectWBCnew(Y= training_sub_rand, coefWBC = meanmeth_subrand)*100

#############################################################################


#calculate R^2 for 72 libraries
cd4rsq72 <- cor(Cell_estimates72[,"CD4T"], truepropscombined[,"CD4T"], method = "pearson")^2
cd8rsq72 <- cor(Cell_estimates72[,"CD8T"], truepropscombined[,"CD8T"], method = "pearson")^2
nkrsq72 <- cor(Cell_estimates72[,"NK"], truepropscombined[,"NK"], method = "pearson")^2
brsq72 <- cor(Cell_estimates72[,"Bcell"], truepropscombined[,"Bcell"], method = "pearson")^2
monorsq72 <- cor(Cell_estimates72[,"Mono"], truepropscombined[,"Monocyte"], method = "pearson")^2
granrsq72 <- cor(Cell_estimates72[,"Gran"], truepropscombined[,"Granulocyte"], method = "pearson")^2
#get the average R squares across all cell types
averagersq72 <- sum(c(cd4rsq72, cd8rsq72, nkrsq72, brsq72, monorsq72, granrsq72))/6

#calculate R^2 for 120 libraries
cd4rsq120 <- cor(Cell_estimates120[,"CD4T"], truepropscombined[,"CD4T"], method = "pearson")^2
cd8rsq120 <- cor(Cell_estimates120[,"CD8T"], truepropscombined[,"CD8T"], method = "pearson")^2
nkrsq120 <- cor(Cell_estimates120[,"NK"], truepropscombined[,"NK"], method = "pearson")^2
brsq120 <- cor(Cell_estimates120[,"Bcell"], truepropscombined[,"Bcell"], method = "pearson")^2
monorsq120 <- cor(Cell_estimates120[,"Mono"], truepropscombined[,"Monocyte"], method = "pearson")^2
granrsq120 <- cor(Cell_estimates120[,"Gran"], truepropscombined[,"Granulocyte"], method = "pearson")^2
#get the average R squares across all cell types
averagersq120 <- sum(c(cd4rsq120, cd8rsq120, nkrsq120, brsq120, monorsq120, granrsq120))/6

#calculate R^2 for 180 libraries
cd4rsq180 <- cor(Cell_estimates180[,"CD4T"], truepropscombined[,"CD4T"], method = "pearson")^2
cd8rsq180 <- cor(Cell_estimates180[,"CD8T"], truepropscombined[,"CD8T"], method = "pearson")^2
nkrsq180 <- cor(Cell_estimates180[,"NK"], truepropscombined[,"NK"], method = "pearson")^2
brsq180 <- cor(Cell_estimates180[,"Bcell"], truepropscombined[,"Bcell"], method = "pearson")^2
monorsq180 <- cor(Cell_estimates180[,"Mono"], truepropscombined[,"Monocyte"], method = "pearson")^2
granrsq180 <- cor(Cell_estimates180[,"Gran"], truepropscombined[,"Granulocyte"], method = "pearson")^2
#get the average R squares across all cell types
averagersq180 <- sum(c(cd4rsq180, cd8rsq180, nkrsq180, brsq180, monorsq180, granrsq180))/6

#calculate R^2 for 240 libraries
cd4rsq240 <- cor(Cell_estimates240[,"CD4T"], truepropscombined[,"CD4T"], method = "pearson")^2
cd8rsq240 <- cor(Cell_estimates240[,"CD8T"], truepropscombined[,"CD8T"], method = "pearson")^2
nkrsq240 <- cor(Cell_estimates240[,"NK"], truepropscombined[,"NK"], method = "pearson")^2
brsq240 <- cor(Cell_estimates240[,"Bcell"], truepropscombined[,"Bcell"], method = "pearson")^2
monorsq240 <- cor(Cell_estimates240[,"Mono"], truepropscombined[,"Monocyte"], method = "pearson")^2
granrsq240 <- cor(Cell_estimates240[,"Gran"], truepropscombined[,"Granulocyte"], method = "pearson")^2
#get the average R squares across all cell types
averagersq240 <- sum(c(cd4rsq240, cd8rsq240, nkrsq240, brsq240, monorsq240, granrsq240))/6


#calculate R^2 for 300 libraries
cd4rsq300 <- cor(Cell_estimates300[,"CD4T"], truepropscombined[,"CD4T"], method = "pearson")^2
cd8rsq300 <- cor(Cell_estimates300[,"CD8T"], truepropscombined[,"CD8T"], method = "pearson")^2
nkrsq300 <- cor(Cell_estimates300[,"NK"], truepropscombined[,"NK"], method = "pearson")^2
brsq300 <- cor(Cell_estimates300[,"Bcell"], truepropscombined[,"Bcell"], method = "pearson")^2
monorsq300 <- cor(Cell_estimates300[,"Mono"], truepropscombined[,"Monocyte"], method = "pearson")^2
granrsq300 <- cor(Cell_estimates300[,"Gran"], truepropscombined[,"Granulocyte"], method = "pearson")^2
#get the average R squares across all cell types
averagersq300 <- sum(c(cd4rsq300, cd8rsq300, nkrsq300, brsq300, monorsq300, granrsq300))/6

#calculate R^2 for 360 libraries
cd4rsq360 <- cor(Cell_estimates360[,"CD4T"], truepropscombined[,"CD4T"], method = "pearson")^2
cd8rsq360 <- cor(Cell_estimates360[,"CD8T"], truepropscombined[,"CD8T"], method = "pearson")^2
nkrsq360 <- cor(Cell_estimates360[,"NK"], truepropscombined[,"NK"], method = "pearson")^2
brsq360 <- cor(Cell_estimates360[,"Bcell"], truepropscombined[,"Bcell"], method = "pearson")^2
monorsq360 <- cor(Cell_estimates360[,"Mono"], truepropscombined[,"Monocyte"], method = "pearson")^2
granrsq360 <- cor(Cell_estimates360[,"Gran"], truepropscombined[,"Granulocyte"], method = "pearson")^2
#get the average R squares across all cell types
averagersq360 <- sum(c(cd4rsq360, cd8rsq360, nkrsq360, brsq360, monorsq360, granrsq360))/6

#calculate R^2 for 180 libraries
cd4rsq540 <- cor(Cell_estimates540[,"CD4T"], truepropscombined[,"CD4T"], method = "pearson")^2
cd8rsq540 <- cor(Cell_estimates540[,"CD8T"], truepropscombined[,"CD8T"], method = "pearson")^2
nkrsq540 <- cor(Cell_estimates540[,"NK"], truepropscombined[,"NK"], method = "pearson")^2
brsq540 <- cor(Cell_estimates540[,"Bcell"], truepropscombined[,"Bcell"], method = "pearson")^2
monorsq540 <- cor(Cell_estimates540[,"Mono"], truepropscombined[,"Monocyte"], method = "pearson")^2
granrsq540 <- cor(Cell_estimates540[,"Gran"], truepropscombined[,"Granulocyte"], method = "pearson")^2
#get the average R squares across all cell types
averagersq540 <- sum(c(cd4rsq540, cd8rsq540, nkrsq540, brsq540, monorsq540, granrsq540))/6

################NOW DO RMSE##############################

#RMSE for 72 library
cd4rmse72 <- rmse(Cell_estimates72[,"CD4T"], truepropscombined[,"CD4T"])
cd8rmse72 <- rmse(Cell_estimates72[,"CD8T"], truepropscombined[,"CD8T"])
nkrmse72 <- rmse(Cell_estimates72[,"NK"], truepropscombined[,"NK"])
brmse72 <- rmse(Cell_estimates72[,"Bcell"], truepropscombined[,"Bcell"])
monormse72 <- rmse(Cell_estimates72[,"Mono"], truepropscombined[,"Monocyte"])
granrmse72 <- rmse(Cell_estimates72[,"Gran"], truepropscombined[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse72 <- sum(c(cd4rmse72, cd8rmse72, nkrmse72, brmse72, monormse72, granrmse72))/6

#RMSE for 120 library
cd4rmse120 <- rmse(Cell_estimates120[,"CD4T"], truepropscombined[,"CD4T"])
cd8rmse120 <- rmse(Cell_estimates120[,"CD8T"], truepropscombined[,"CD8T"])
nkrmse120 <- rmse(Cell_estimates120[,"NK"], truepropscombined[,"NK"])
brmse120 <- rmse(Cell_estimates120[,"Bcell"], truepropscombined[,"Bcell"])
monormse120 <- rmse(Cell_estimates120[,"Mono"], truepropscombined[,"Monocyte"])
granrmse120 <- rmse(Cell_estimates120[,"Gran"], truepropscombined[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse120 <- sum(c(cd4rmse120, cd8rmse120, nkrmse120, brmse120, monormse120, granrmse120))/6

#RMSE for 180 library
cd4rmse180 <- rmse(Cell_estimates180[,"CD4T"], truepropscombined[,"CD4T"])
cd8rmse180 <- rmse(Cell_estimates180[,"CD8T"], truepropscombined[,"CD8T"])
nkrmse180 <- rmse(Cell_estimates180[,"NK"], truepropscombined[,"NK"])
brmse180 <- rmse(Cell_estimates180[,"Bcell"], truepropscombined[,"Bcell"])
monormse180 <- rmse(Cell_estimates180[,"Mono"], truepropscombined[,"Monocyte"])
granrmse180 <- rmse(Cell_estimates180[,"Gran"], truepropscombined[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse180 <- sum(c(cd4rmse180, cd8rmse180, nkrmse180, brmse180, monormse180, granrmse180))/6

#RMSE for 240 library
cd4rmse240 <- rmse(Cell_estimates240[,"CD4T"], truepropscombined[,"CD4T"])
cd8rmse240 <- rmse(Cell_estimates240[,"CD8T"], truepropscombined[,"CD8T"])
nkrmse240 <- rmse(Cell_estimates240[,"NK"], truepropscombined[,"NK"])
brmse240 <- rmse(Cell_estimates240[,"Bcell"], truepropscombined[,"Bcell"])
monormse240 <- rmse(Cell_estimates240[,"Mono"], truepropscombined[,"Monocyte"])
granrmse240 <- rmse(Cell_estimates240[,"Gran"], truepropscombined[,6])
#get the average RMSE across all cell types
averagermse240 <- sum(c(cd4rmse240, cd8rmse240, nkrmse240, brmse240, monormse240, granrmse240))/6

#RMSE for 300 library
cd4rmse300 <- rmse(Cell_estimates300[,"CD4T"], truepropscombined[,"CD4T"])
cd8rmse300 <- rmse(Cell_estimates300[,"CD8T"], truepropscombined[,"CD8T"])
nkrmse300 <- rmse(Cell_estimates300[,"NK"], truepropscombined[,"NK"])
brmse300 <- rmse(Cell_estimates300[,"Bcell"], truepropscombined[,"Bcell"])
monormse300 <- rmse(Cell_estimates300[,"Mono"], truepropscombined[,"Monocyte"])
granrmse300 <- rmse(Cell_estimates300[,"Gran"], truepropscombined[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse300 <- sum(c(cd4rmse300, cd8rmse300, nkrmse300, brmse300, monormse300, granrmse300))/6

#RMSE for 360 library
cd4rmse360 <- rmse(Cell_estimates360[,"CD4T"], truepropscombined[,"CD4T"])
cd8rmse360 <- rmse(Cell_estimates360[,"CD8T"], truepropscombined[,"CD8T"])
nkrmse360 <- rmse(Cell_estimates360[,"NK"], truepropscombined[,"NK"])
brmse360 <- rmse(Cell_estimates360[,"Bcell"], truepropscombined[,"Bcell"])
monormse360 <- rmse(Cell_estimates360[,"Mono"], truepropscombined[,"Monocyte"])
granrmse360 <- rmse(Cell_estimates360[,"Gran"], truepropscombined[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse360 <- sum(c(cd4rmse360, cd8rmse360, nkrmse360, brmse360, monormse360, granrmse360))/6

#RMSE for 540 library
cd4rmse540 <- rmse(Cell_estimates540[,"CD4T"], truepropscombined[,"CD4T"])
cd8rmse540 <- rmse(Cell_estimates540[,"CD8T"], truepropscombined[,"CD8T"])
nkrmse540 <- rmse(Cell_estimates540[,"NK"], truepropscombined[,"NK"])
brmse540 <- rmse(Cell_estimates540[,"Bcell"], truepropscombined[,"Bcell"])
monormse540 <- rmse(Cell_estimates540[,"Mono"], truepropscombined[,"Monocyte"])
granrmse540 <- rmse(Cell_estimates540[,"Gran"], truepropscombined[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse540 <- sum(c(cd4rmse540, cd8rmse540, nkrmse540, brmse540, monormse540, granrmse540))/6



#get info for random library

cd4rsqrand <- cor(Cell_estimatesrand[,"CD4T"], truepropscombined[,"CD4T"], method = "pearson")^2
cd8rsqrand <- cor(Cell_estimatesrand[,"CD8T"], truepropscombined[,"CD8T"], method = "pearson")^2
nkrsqrand <- cor(Cell_estimatesrand[,"NK"], truepropscombined[,"NK"], method = "pearson")^2
brsqrand <- cor(Cell_estimatesrand[,"Bcell"], truepropscombined[,"Bcell"], method = "pearson")^2
monorsqrand <- cor(Cell_estimatesrand[,"Mono"], truepropscombined[,"Monocyte"], method = "pearson")^2
granrsqrand <- cor(Cell_estimatesrand[,"Gran"], truepropscombined[,"Granulocyte"], method = "pearson")^2
#get the average R squares across all cell types
averagersqrand <- sum(c(cd4rsqrand, cd8rsqrand, nkrsqrand, brsqrand, monorsqrand, granrsqrand))/6

cd4rmserand <- rmse(Cell_estimatesrand[,"CD4T"], truepropscombined[,"CD4T"])
cd8rmserand <- rmse(Cell_estimatesrand[,"CD8T"], truepropscombined[,"CD8T"])
nkrmserand <- rmse(Cell_estimatesrand[,"NK"], truepropscombined[,"NK"])
brmserand <- rmse(Cell_estimatesrand[,"Bcell"], truepropscombined[,"Bcell"])
monormserand <- rmse(Cell_estimatesrand[,"Mono"], truepropscombined[,"Monocyte"])
granrmserand <- rmse(Cell_estimatesrand[,"Gran"], truepropscombined[,"Granulocyte"])
#get the average RMSE across all cell types
averagermserand <- sum(c(cd4rmserand, cd8rmserand, nkrmserand, brmserand, monormserand, granrmserand))/6



########Now get information for t tests#################


#get library using t tests
candobject1 = CandidateDMRFinder.v2(uniquecelltypes, referenceBetas, referencePd, M = 6)
candobject2 = CandidateDMRFinder.v2(uniquecelltypes, referenceBetas, referencePd, M = 10)
candobject3 = CandidateDMRFinder.v2(uniquecelltypes, referenceBetas, referencePd, M = 15)
candobject4 = CandidateDMRFinder.v2(uniquecelltypes, referenceBetas, referencePd, M = 20)
candobject5 = CandidateDMRFinder.v2(uniquecelltypes, referenceBetas, referencePd, M = 25)
candobject6 = CandidateDMRFinder.v2(uniquecelltypes, referenceBetas, referencePd, M = 30)
candobject7 = CandidateDMRFinder.v2(uniquecelltypes, referenceBetas, referencePd, M = 45)


tstats72 <- candobject1$candidateSet
tstats120 <- candobject2$candidateSet
tstats180 <- candobject3$candidateSet
tstats240 <- candobject4$candidateSet
tstats300 <- candobject5$candidateSet
tstats360 <- candobject6$candidateSet
tstats540 <- candobject7$candidateSet

#save the CpGs
save(tstats120, file="Legacy120.RData")
save(tstats72, file="Legacy72.RData")
save(tstats180, file="Legacy180.RData")
save(tstats240, file="Legacy240.RData")
save(tstats300, file="Legacy300.RData")
save(tstats360, file="Legacy360.RData")
save(tstats540, file="Legacy540.RData")


#use full training data set and subset
training_sub72_t <- combinedbetas[tstats72,]
training_sub120_t <- combinedbetas[tstats120,]
training_sub180_t <- combinedbetas[tstats180,]
training_sub240_t <- combinedbetas[tstats240,]
training_sub300_t <- combinedbetas[tstats300,]
training_sub360_t <- combinedbetas[tstats360,]
training_sub540_t <- combinedbetas[tstats540,]

#get the mean meth values
meanmeth_sub72_t <- meanmeth[tstats72,]
meanmeth_sub120_t <- meanmeth[tstats120,]
meanmeth_sub180_t <- meanmeth[tstats180,]
meanmeth_sub240_t <- meanmeth[tstats240,]
meanmeth_sub300_t <- meanmeth[tstats300,]
meanmeth_sub360_t <- meanmeth[tstats360,]
meanmeth_sub540_t <- meanmeth[tstats540,]

#get predicted cell proportions for each cell types using training data
Cell_estimates72_t <- projectWBCnew(Y= training_sub72_t, coefWBC = meanmeth_sub72_t)*100
Cell_estimates120_t <- projectWBCnew(Y= training_sub120_t, coefWBC = meanmeth_sub120_t)*100
Cell_estimates180_t <- projectWBCnew(Y= training_sub180_t, coefWBC = meanmeth_sub180_t)*100
Cell_estimates240_t <- projectWBCnew(Y= training_sub240_t, coefWBC = meanmeth_sub240_t)*100
Cell_estimates300_t <- projectWBCnew(Y= training_sub300_t, coefWBC = meanmeth_sub300_t)*100
Cell_estimates360_t <- projectWBCnew(Y= training_sub360_t, coefWBC = meanmeth_sub360_t)*100
Cell_estimates540_t <- projectWBCnew(Y= training_sub540_t, coefWBC = meanmeth_sub540_t)*100

#calculate R^2 for 72 libraries
cd4rsq72_t <- cor(Cell_estimates72_t[,"CD4T"], truepropscombined[,"CD4T"], method = "pearson")^2
cd8rsq72_t <- cor(Cell_estimates72_t[,"CD8T"], truepropscombined[,"CD8T"], method = "pearson")^2
nkrsq72_t <- cor(Cell_estimates72_t[,"NK"], truepropscombined[,"NK"], method = "pearson")^2
brsq72_t <- cor(Cell_estimates72_t[,"Bcell"], truepropscombined[,"Bcell"], method = "pearson")^2
monorsq72_t <- cor(Cell_estimates72_t[,"Mono"], truepropscombined[,"Monocyte"], method = "pearson")^2
granrsq72_t <- cor(Cell_estimates72_t[,"Gran"], truepropscombined[,"Granulocyte"], method = "pearson")^2
#get the average R squares across all cell types
averagersq72_t <- sum(c(cd4rsq72_t, cd8rsq72_t, nkrsq72_t, brsq72_t, monorsq72_t, granrsq72_t))/6

#calculate R^2 for 120 libraries
cd4rsq120_t <- cor(Cell_estimates120_t[,"CD4T"], truepropscombined[,"CD4T"], method = "pearson")^2
cd8rsq120_t <- cor(Cell_estimates120_t[,"CD8T"], truepropscombined[,"CD8T"], method = "pearson")^2
nkrsq120_t <- cor(Cell_estimates120_t[,"NK"], truepropscombined[,"NK"], method = "pearson")^2
brsq120_t <- cor(Cell_estimates120_t[,"Bcell"], truepropscombined[,"Bcell"], method = "pearson")^2
monorsq120_t <- cor(Cell_estimates120_t[,"Mono"], truepropscombined[,"Monocyte"], method = "pearson")^2
granrsq120_t <- cor(Cell_estimates120_t[,"Gran"], truepropscombined[,"Granulocyte"], method = "pearson")^2
#get the average R squares across all cell types
averagersq120_t <- sum(c(cd4rsq120_t, cd8rsq120_t, nkrsq120_t, brsq120_t, monorsq120_t, granrsq120_t))/6

#calculate R^2 for 180 libraries
cd4rsq180_t <- cor(Cell_estimates180_t[,"CD4T"], truepropscombined[,"CD4T"], method = "pearson")^2
cd8rsq180_t <- cor(Cell_estimates180_t[,"CD8T"], truepropscombined[,"CD8T"], method = "pearson")^2
nkrsq180_t <- cor(Cell_estimates180_t[,"NK"], truepropscombined[,"NK"], method = "pearson")^2
brsq180_t <- cor(Cell_estimates180_t[,"Bcell"], truepropscombined[,"Bcell"], method = "pearson")^2
monorsq180_t <- cor(Cell_estimates180_t[,"Mono"], truepropscombined[,"Monocyte"], method = "pearson")^2
granrsq180_t <- cor(Cell_estimates180_t[,"Gran"], truepropscombined[,"Granulocyte"], method = "pearson")^2
#get the average R squares across all cell types
averagersq180_t <- sum(c(cd4rsq180_t, cd8rsq180_t, nkrsq180_t, brsq180_t, monorsq180_t, granrsq180_t))/6

#calculate R^2 for 240 libraries
cd4rsq240_t <- cor(Cell_estimates240_t[,"CD4T"], truepropscombined[,"CD4T"], method = "pearson")^2
cd8rsq240_t <- cor(Cell_estimates240_t[,"CD8T"], truepropscombined[,"CD8T"], method = "pearson")^2
nkrsq240_t <- cor(Cell_estimates240_t[,"NK"], truepropscombined[,"NK"], method = "pearson")^2
brsq240_t <- cor(Cell_estimates240_t[,"Bcell"], truepropscombined[,"Bcell"], method = "pearson")^2
monorsq240_t <- cor(Cell_estimates240_t[,"Mono"], truepropscombined[,"Monocyte"], method = "pearson")^2
granrsq240_t <- cor(Cell_estimates240_t[,"Gran"], truepropscombined[,"Granulocyte"], method = "pearson")^2
#get the average R squares across all cell types
averagersq240_t <- sum(c(cd4rsq240_t, cd8rsq240_t, nkrsq240_t, brsq240_t, monorsq240_t, granrsq240_t))/6

#calculate R^2 for 300 libraries
cd4rsq300_t <- cor(Cell_estimates300_t[,"CD4T"], truepropscombined[,"CD4T"], method = "pearson")^2
cd8rsq300_t <- cor(Cell_estimates300_t[,"CD8T"], truepropscombined[,"CD8T"], method = "pearson")^2
nkrsq300_t <- cor(Cell_estimates300_t[,"NK"], truepropscombined[,"NK"], method = "pearson")^2
brsq300_t <- cor(Cell_estimates300_t[,"Bcell"], truepropscombined[,"Bcell"], method = "pearson")^2
monorsq300_t <- cor(Cell_estimates300_t[,"Mono"], truepropscombined[,"Monocyte"], method = "pearson")^2
granrsq300_t <- cor(Cell_estimates300_t[,"Gran"], truepropscombined[,"Granulocyte"], method = "pearson")^2
#get the average R squares across all cell types
averagersq300_t <- sum(c(cd4rsq300_t, cd8rsq300_t, nkrsq300_t, brsq300_t, monorsq300_t, granrsq300_t))/6

#calculate R^2 for 360 libraries
cd4rsq360_t <- cor(Cell_estimates360_t[,"CD4T"], truepropscombined[,"CD4T"], method = "pearson")^2
cd8rsq360_t <- cor(Cell_estimates360_t[,"CD8T"], truepropscombined[,"CD8T"], method = "pearson")^2
nkrsq360_t <- cor(Cell_estimates360_t[,"NK"], truepropscombined[,"NK"], method = "pearson")^2
brsq360_t <- cor(Cell_estimates360_t[,"Bcell"], truepropscombined[,"Bcell"], method = "pearson")^2
monorsq360_t <- cor(Cell_estimates360_t[,"Mono"], truepropscombined[,"Monocyte"], method = "pearson")^2
granrsq360_t <- cor(Cell_estimates360_t[,"Gran"], truepropscombined[,"Granulocyte"], method = "pearson")^2
#get the average R squares across all cell types
averagersq360_t <- sum(c(cd4rsq360_t, cd8rsq360_t, nkrsq360_t, brsq360_t, monorsq360_t, granrsq360_t))/6

#calculate R^2 for 540 libraries
cd4rsq540_t <- cor(Cell_estimates540_t[,"CD4T"], truepropscombined[,"CD4T"], method = "pearson")^2
cd8rsq540_t <- cor(Cell_estimates540_t[,"CD8T"], truepropscombined[,"CD8T"], method = "pearson")^2
nkrsq540_t <- cor(Cell_estimates540_t[,"NK"], truepropscombined[,"NK"], method = "pearson")^2
brsq540_t <- cor(Cell_estimates540_t[,"Bcell"], truepropscombined[,"Bcell"], method = "pearson")^2
monorsq540_t <- cor(Cell_estimates540_t[,"Mono"], truepropscombined[,"Monocyte"], method = "pearson")^2
granrsq540_t <- cor(Cell_estimates540_t[,"Gran"], truepropscombined[,"Granulocyte"], method = "pearson")^2
#get the average R squares across all cell types
averagersq540_t <- sum(c(cd4rsq540_t, cd8rsq540_t, nkrsq540_t, brsq540_t, monorsq540_t, granrsq540_t))/6



#now get rmse for 72 library
cd4rmse72_t <- rmse(Cell_estimates72_t[,"CD4T"], truepropscombined[,"CD4T"])
cd8rmse72_t <- rmse(Cell_estimates72_t[,"CD8T"], truepropscombined[,"CD8T"])
nkrmse72_t <- rmse(Cell_estimates72_t[,"NK"], truepropscombined[,"NK"])
brmse72_t <- rmse(Cell_estimates72_t[,"Bcell"], truepropscombined[,"Bcell"])
monormse72_t <- rmse(Cell_estimates72_t[,"Mono"], truepropscombined[,"Monocyte"])
granrmse72_t <- rmse(Cell_estimates72_t[,"Gran"], truepropscombined[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse72_t <- sum(c(cd4rmse72_t, cd8rmse72_t, nkrmse72_t, brmse72_t, monormse72_t, granrmse72_t))/6

#now get rmse for 120 library
cd4rmse120_t <- rmse(Cell_estimates120_t[,"CD4T"], truepropscombined[,"CD4T"])
cd8rmse120_t <- rmse(Cell_estimates120_t[,"CD8T"], truepropscombined[,"CD8T"])
nkrmse120_t <- rmse(Cell_estimates120_t[,"NK"], truepropscombined[,"NK"])
brmse120_t <- rmse(Cell_estimates120_t[,"Bcell"], truepropscombined[,"Bcell"])
monormse120_t <- rmse(Cell_estimates120_t[,"Mono"], truepropscombined[,"Monocyte"])
granrmse120_t <- rmse(Cell_estimates120_t[,"Gran"], truepropscombined[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse120_t <- sum(c(cd4rmse120_t, cd8rmse120_t, nkrmse120_t, brmse120_t, monormse120_t, granrmse120_t))/6

#now get rmse for 180 library
cd4rmse180_t <- rmse(Cell_estimates180_t[,"CD4T"], truepropscombined[,"CD4T"])
cd8rmse180_t <- rmse(Cell_estimates180_t[,"CD8T"], truepropscombined[,"CD8T"])
nkrmse180_t <- rmse(Cell_estimates180_t[,"NK"], truepropscombined[,"NK"])
brmse180_t <- rmse(Cell_estimates180_t[,"Bcell"], truepropscombined[,"Bcell"])
monormse180_t <- rmse(Cell_estimates180_t[,"Mono"], truepropscombined[,"Monocyte"])
granrmse180_t <- rmse(Cell_estimates180_t[,"Gran"], truepropscombined[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse180_t <- sum(c(cd4rmse180_t, cd8rmse180_t, nkrmse180_t, brmse180_t, monormse180_t, granrmse180_t))/6

#now get rmse for 240 library
cd4rmse240_t <- rmse(Cell_estimates240_t[,"CD4T"], truepropscombined[,"CD4T"])
cd8rmse240_t <- rmse(Cell_estimates240_t[,"CD8T"], truepropscombined[,"CD8T"])
nkrmse240_t <- rmse(Cell_estimates240_t[,"NK"], truepropscombined[,"NK"])
brmse240_t <- rmse(Cell_estimates240_t[,"Bcell"], truepropscombined[,"Bcell"])
monormse240_t <- rmse(Cell_estimates240_t[,"Mono"], truepropscombined[,"Monocyte"])
granrmse240_t <- rmse(Cell_estimates240_t[,"Gran"], truepropscombined[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse240_t <- sum(c(cd4rmse240_t, cd8rmse240_t, nkrmse240_t, brmse240_t, monormse240_t, granrmse240_t))/6

#now get rmse for 300 library
cd4rmse300_t <- rmse(Cell_estimates300_t[,"CD4T"], truepropscombined[,"CD4T"])
cd8rmse300_t <- rmse(Cell_estimates300_t[,"CD8T"], truepropscombined[,"CD8T"])
nkrmse300_t <- rmse(Cell_estimates300_t[,"NK"], truepropscombined[,"NK"])
brmse300_t <- rmse(Cell_estimates300_t[,"Bcell"], truepropscombined[,"Bcell"])
monormse300_t <- rmse(Cell_estimates300_t[,"Mono"], truepropscombined[,"Monocyte"])
granrmse300_t <- rmse(Cell_estimates300_t[,"Gran"], truepropscombined[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse300_t <- sum(c(cd4rmse300_t, cd8rmse300_t, nkrmse300_t, brmse300_t, monormse300_t, granrmse300_t))/6

#now get rmse for 360 library
cd4rmse360_t <- rmse(Cell_estimates360_t[,"CD4T"], truepropscombined[,"CD4T"])
cd8rmse360_t <- rmse(Cell_estimates360_t[,"CD8T"], truepropscombined[,"CD8T"])
nkrmse360_t <- rmse(Cell_estimates360_t[,"NK"], truepropscombined[,"NK"])
brmse360_t <- rmse(Cell_estimates360_t[,"Bcell"], truepropscombined[,"Bcell"])
monormse360_t <- rmse(Cell_estimates360_t[,"Mono"], truepropscombined[,"Monocyte"])
granrmse360_t <- rmse(Cell_estimates360_t[,"Gran"], truepropscombined[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse360_t <- sum(c(cd4rmse360_t, cd8rmse360_t, nkrmse360_t, brmse360_t, monormse360_t, granrmse360_t))/6

#now get rmse for 540 library
cd4rmse540_t <- rmse(Cell_estimates540_t[,"CD4T"], truepropscombined[,"CD4T"])
cd8rmse540_t <- rmse(Cell_estimates540_t[,"CD8T"], truepropscombined[,"CD8T"])
nkrmse540_t <- rmse(Cell_estimates540_t[,"NK"], truepropscombined[,"NK"])
brmse540_t <- rmse(Cell_estimates540_t[,"Bcell"], truepropscombined[,"Bcell"])
monormse540_t <- rmse(Cell_estimates540_t[,"Mono"], truepropscombined[,"Monocyte"])
granrmse540_t <- rmse(Cell_estimates540_t[,"Gran"], truepropscombined[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse540_t <- sum(c(cd4rmse540_t, cd8rmse540_t, nkrmse540_t, brmse540_t, monormse540_t, granrmse540_t))/6


########NOW PUT INFO INTO TABLES#############

# R^2

#Average across all 6 cell types
newmetavg <- c(averagersq72, averagersq120, averagersq180,averagersq240,averagersq300, averagersq360,
               averagersq540)
tstatavgrsq <- c(averagersq72_t, averagersq120_t,averagersq180_t,averagersq240_t,
              averagersq300_t,averagersq360_t,averagersq540_t)

compareavgrsq <- cbind(newmetavg, tstatavgrsq)
rownames(compareavgrsq) <- c("72 Cpgs", "120 Cpgs","180 Cpgs","240 Cpgs","300 Cpgs","360 Cpgs","540 Cpgs")

#rmse 

#Average across all 6 tell types
newmetavgrmse <- c(averagermse72, averagermse120, averagermse180,averagermse240,averagermse300, averagermse360,
                   averagermse540)
tstatavgrmse <- c(averagermse72_t, averagermse120_t,averagermse180_t,averagermse240_t,
                  averagermse300_t,averagermse360_t,averagermse540_t)

compareavgrmse <- cbind(newmetavgrmse, tstatavgrmse)
rownames(compareavgrmse) <- c("72 Cpgs", "120 Cpgs","180 Cpgs","240 Cpgs","300 Cpgs","360 Cpgs","540 Cpgs")

#results for random library
randomlib_results <- c(averagersqrand, averagermserand)
names(randomlib_results) <- c("R^2", "RMSE")


#forcd4
newmetcd4rsq <- c(cd4rsq72, cd4rsq120,cd4rsq180,cd4rsq240,cd4rsq300,cd4rsq360,cd4rsq540)
tstatcd4rsq <- c(cd4rsq72_t,cd4rsq120_t,cd4rsq180_t,cd4rsq240_t,cd4rsq300_t,cd4rsq360_t,cd4rsq540_t)
newmetcd4rmse <- c(cd4rmse72, cd4rmse120,cd4rmse180,cd4rmse240,cd4rmse300,cd4rmse360,cd4rmse540)
tstatcd4rmse <- c(cd4rmse72_t,cd4rmse120_t,cd4rmse180_t,cd4rmse240_t,cd4rmse300_t,cd4rmse360_t,cd4rmse540_t)

comparecd4 <- cbind(newmetcd4rsq, tstatcd4rsq, newmetcd4rmse,tstatcd4rmse)
rownames(comparecd4) <- c("72 Cpgs", "120 Cpgs","180 Cpgs","240 Cpgs","300 Cpgs","360 Cpgs","540 Cpgs")

#for cd8
newmetcd8rsq <- c(cd8rsq72, cd8rsq120,cd8rsq180,cd8rsq240,cd8rsq300,cd8rsq360,cd8rsq540)
tstatcd8rsq <- c(cd8rsq72_t,cd8rsq120_t,cd8rsq180_t,cd8rsq240_t,cd8rsq300_t,cd8rsq360_t,cd8rsq540_t)
newmetcd8rmse <- c(cd8rmse72, cd8rmse120,cd8rmse180,cd8rmse240,cd8rmse300,cd8rmse360,cd8rmse540)
tstatcd8rmse <- c(cd8rmse72_t,cd8rmse120_t,cd8rmse180_t,cd8rmse240_t,cd8rmse300_t,cd8rmse360_t,cd8rmse540_t)

comparecd8 <- cbind(newmetcd8rsq, tstatcd8rsq, newmetcd8rmse,tstatcd8rmse)
rownames(comparecd8) <- c("72 Cpgs", "120 Cpgs","180 Cpgs","240 Cpgs","300 Cpgs","360 Cpgs","540 Cpgs")

#for nk cells
newmetnkrsq <- c(nkrsq72, nkrsq120,nkrsq180,nkrsq240,nkrsq300,nkrsq360,nkrsq540)
tstatnkrsq <- c(nkrsq72_t,nkrsq120_t,nkrsq180_t,nkrsq240_t,nkrsq300_t,nkrsq360_t,nkrsq540_t)
newmetnkrmse <- c(nkrmse72, nkrmse120,nkrmse180,nkrmse240,nkrmse300,nkrmse360,nkrmse540)
tstatnkrmse <- c(nkrmse72_t,nkrmse120_t,nkrmse180_t,nkrmse240_t,nkrmse300_t,nkrmse360_t,nkrmse540_t)

comparenk <- cbind(newmetnkrsq, tstatnkrsq, newmetnkrmse,tstatnkrmse)
rownames(comparenk) <- c("72 Cpgs", "120 Cpgs","180 Cpgs","240 Cpgs","300 Cpgs","360 Cpgs","540 Cpgs")

#for b cells
newmetbrsq <- c(brsq72, brsq120,brsq180,brsq240,brsq300,brsq360,brsq540)
tstatbrsq <- c(brsq72_t,brsq120_t,brsq180_t,brsq240_t,brsq300_t,brsq360_t,brsq540_t)
newmetbrmse <- c(brmse72, brmse120,brmse180,brmse240,brmse300,brmse360,brmse540)
tstatbrmse <- c(brmse72_t,brmse120_t,brmse180_t,brmse240_t,brmse300_t,brmse360_t,brmse540_t)

compareb <- cbind(newmetbrsq, tstatbrsq, newmetbrmse,tstatbrmse)
rownames(compareb) <- c("72 Cpgs", "120 Cpgs","180 Cpgs","240 Cpgs","300 Cpgs","360 Cpgs","540 Cpgs")

#for monocytes
newmetmonorsq <- c(monorsq72, monorsq120,monorsq180,monorsq240,monorsq300,monorsq360,monorsq540)
tstatmonorsq <- c(monorsq72_t,monorsq120_t,monorsq180_t,monorsq240_t,monorsq300_t,monorsq360_t,monorsq540_t)
newmetmonormse <- c(monormse72, monormse120,monormse180,monormse240,monormse300,monormse360,monormse540)
tstatmonormse <- c(monormse72_t,monormse120_t,monormse180_t,monormse240_t,monormse300_t,monormse360_t,monormse540_t)

comparemono <- cbind(newmetmonorsq, tstatmonorsq, newmetmonormse,tstatmonormse)
rownames(comparemono) <- c("72 Cpgs", "120 Cpgs","180 Cpgs","240 Cpgs","300 Cpgs","360 Cpgs","540 Cpgs")

#for granulacytes
newmetgranrsq <- c(granrsq72, granrsq120,granrsq180,granrsq240,granrsq300,granrsq360,granrsq540)
tstatgranrsq <- c(granrsq72_t,granrsq120_t,granrsq180_t,granrsq240_t,granrsq300_t,granrsq360_t,granrsq540_t)
newmetgranrmse <- c(granrmse72, granrmse120,granrmse180,granrmse240,granrmse300,granrmse360,granrmse540)
tstatgranrmse <- c(granrmse72_t,granrmse120_t,granrmse180_t,granrmse240_t,granrmse300_t,granrmse360_t,granrmse540_t)

comparegran <- cbind(newmetgranrsq, tstatgranrsq, newmetgranrmse,tstatgranrmse)
rownames(comparegran) <- c("72 Cpgs", "120 Cpgs","180 Cpgs","240 Cpgs","300 Cpgs","360 Cpgs","540 Cpgs")

#random library results by cell type
randrsq <- c(cd4rsqrand, cd8rsqrand, brsqrand, nkrsqrand, monorsqrand, granrsqrand)
randrmse <- c(cd4rmserand, cd8rmserand, brmserand, nkrmserand, monormserand, granrmserand)
comparerand <- cbind(randrsq, randrmse)
rownames(comparerand) <- c("CD4", "CD8", "B", "NK", "Mono", "Gran")
colnames(comparerand) <- c("R^2", "RMSE")

#############HEAT MAPS#################

library(gplots)

col.name.order <- c(13, 14, 15, 31, 32, 33, 1, 2, 3, 19, 20, 21, 10, 11, 12, 28, 29, 30, 16, 17, 18,
                    34, 35, 36, 4, 5, 6, 22, 23, 24, 7, 8 , 9, 25, 26, 27)



#new metric heat maps
rand_mat_72 <- RefBetas[library72, col.name.order]
rand_mat_120 <- RefBetas[library120, col.name.order]
rand_mat_180 <- RefBetas[library180, col.name.order]
rand_mat_240 <- RefBetas[library240, col.name.order]
rand_mat_300 <- RefBetas[library300, col.name.order]
rand_mat_360 <- RefBetas[library360, col.name.order]
rand_mat_540 <- RefBetas[library540, col.name.order]

m <- c(5,5)
heatmap.2(rand_mat_72, col = colorRampPalette(c("yellow", "black", "blue"))(32), 
          trace = "n", dendrogram = "row", Rowv = T, Colv = F,margins = m, main= "New Metric 72 Library")

heatmap.2(rand_mat_120, col = colorRampPalette(c("yellow", "black", "blue"))(32), 
          trace = "n", dendrogram = "row", Rowv = T, Colv = F,margins = m, main= "New Metric 120 Library")

heatmap.2(rand_mat_180, col = colorRampPalette(c("yellow", "black", "blue"))(32), 
          trace = "n", dendrogram = "row", Rowv = T, Colv = F,margins = m, main= "New Metric 180 Library")

heatmap.2(rand_mat_240, col = colorRampPalette(c("yellow", "black", "blue"))(32), 
          trace = "n", dendrogram = "row", Rowv = T, Colv = F,margins = m, main= "New Metric 240 Library")

heatmap.2(rand_mat_300, col = colorRampPalette(c("yellow", "black", "blue"))(32), 
          trace = "n", dendrogram = "row", Rowv = T, Colv = F,margins = m, main= "New Metric 300 Library")

heatmap.2(rand_mat_360, col = colorRampPalette(c("yellow", "black", "blue"))(32), 
          trace = "n", dendrogram = "row", Rowv = T, Colv = F,margins = m, main= "New Metric 360 Library")

heatmap.2(rand_mat_540, col = colorRampPalette(c("yellow", "black", "blue"))(32), 
          trace = "n", dendrogram = "row", Rowv = T, Colv = F,margins = m, main= "New Metric 540 Library")

#t stat heat maps

t_mat_72 <- RefBetas[tstats72, col.name.order]
t_mat_120 <- RefBetas[tstats120, col.name.order]
t_mat_180 <- RefBetas[tstats180, col.name.order]
t_mat_240 <- RefBetas[tstats240, col.name.order]
t_mat_300 <- RefBetas[tstats300, col.name.order]
t_mat_360 <- RefBetas[tstats360, col.name.order]
t_mat_540 <- RefBetas[tstats540, col.name.order]


heatmap.2(t_mat_72, col = colorRampPalette(c("yellow", "black", "blue"))(32), 
          trace = "n", dendrogram = "row", Rowv = T, Colv = F,margins = m, main= "T Stat 72 Library")

heatmap.2(t_mat_120, col = colorRampPalette(c("yellow", "black", "blue"))(32), 
          trace = "n", dendrogram = "row", Rowv = T, Colv = F,margins = m, main= "T Stat 120 Library")

heatmap.2(t_mat_180, col = colorRampPalette(c("yellow", "black", "blue"))(32), 
          trace = "n", dendrogram = "row", Rowv = T, Colv = F,margins = m, main= "T Stat 180 Library")

heatmap.2(t_mat_240, col = colorRampPalette(c("yellow", "black", "blue"))(32), 
          trace = "n", dendrogram = "row", Rowv = T, Colv = F,margins = m, main= "T Stat 240 Library")

heatmap.2(t_mat_300, col = colorRampPalette(c("yellow", "black", "blue"))(32), 
          trace = "n", dendrogram = "row", Rowv = T, Colv = F,margins = m, main= "T Stat 300 Library")

heatmap.2(t_mat_360, col = colorRampPalette(c("yellow", "black", "blue"))(32), 
          trace = "n", dendrogram = "row", Rowv = T, Colv = F,margins = m, main= "T Stat 360 Library")

heatmap.2(t_mat_540, col = colorRampPalette(c("yellow", "black", "blue"))(32), 
          trace = "n", dendrogram = "row", Rowv = T, Colv = F,margins = m, main= "T Stat 540 Library")



##################################################################################
# PREDICTED VS ACTUAL PLOTS

################################################################################
#fixed library size of 300

library(ggplot2)
library(ggrepel)
library(gridExtra)

#plot for CD4
cd4pred <- data.frame(Cell_estimates300[,"CD4T"], truepropscombined[,"CD4T"])
lbl <- paste("R^2 == ", 0.906)
cd4plot <- ggplot(data = cd4pred, aes(x=truepropscombined[,"CD4T"], y=Cell_estimates300[,"CD4T"])) + 
  geom_point(color="blue", size=2) + 
  ggtitle("CD4 Cells") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                        axis.title.x = element_text(size=14, face="bold"), 
                         axis.title.y = element_text(size=14, face="bold"),
                        axis.text.x = element_text(face="bold",size=12),
                                                                  axis.text.y = element_text(face="bold",size=12)) +
  xlab("True Proportions") + ylab("Predicted Proportions") + 
  geom_abline(size=1) + 
  annotate("text",x=12, y=24,label=lbl, parse=TRUE) + 
  annotate("text",x=12, y=22,label="RMSE = 1.90")



#plot for CD8

cd8pred <- data.frame(Cell_estimates300[,"CD8T"], truepropscombined[,"CD8T"])
lbl <- paste("R^2 == ", 0.953)
cd8plot <- ggplot(data = cd8pred, aes(x=truepropscombined[,"CD8T"], y=Cell_estimates300[,"CD8T"])) + 
  geom_point(color="red", size=2) + ggtitle("CD8 Cells") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold"),
        axis.title.x = element_text(size=14, face="bold"), 
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12)) +
  xlab("True Proportions") + ylab("Predicted Proportions") + 
  geom_abline(size=1) + 
  annotate("text",x=11, y=27,label=lbl, parse=TRUE) + 
  annotate("text",x=11, y=24,label="RMSE = 3.20")



#plot for B cell
bcellpred <- data.frame(Cell_estimates300[,"Bcell"], truepropscombined[,"Bcell"])
lbl <- paste("R^2 == ", 0.993)
bcellplot <- ggplot(data = bcellpred, aes(x=truepropscombined[,"Bcell"], y=Cell_estimates300[,"Bcell"])) + 
  geom_point(color="dark green", size=2) + ggtitle("B Cells") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold"),
        axis.title.x = element_text(size=14, face="bold"), 
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12)) +
  xlab("True Proportions") + ylab("Predicted Proportions") + 
  geom_abline(size=1) + 
  annotate("text",x=9, y=22,label=lbl, parse=TRUE) + 
  annotate("text",x=9, y=20,label="RMSE = 0.624")



#plot for NK
nkpred <- data.frame(Cell_estimates300[,"NK"], truepropscombined[,"NK"])
lbl <- paste("R^2 == ", 0.820)
nkplot <- ggplot(data = nkpred, aes(x=truepropscombined[,"NK"], y=Cell_estimates300[,"NK"])) + 
  geom_point(color="purple", size=2) + ggtitle("Natural Killer Cells") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold"),
        axis.title.x = element_text(size=14, face="bold"), 
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12)) +
  xlab("True Proportions") + ylab("Predicted Proportions") + 
  geom_abline(size=1) + 
  annotate("text",x=7, y=23,label=lbl, parse=TRUE) + 
  annotate("text",x=7, y=21,label="RMSE = 3.49")



#plot for Monocytes
monopred <- data.frame(Cell_estimates300[,"Mono"], truepropscombined[,"Monocyte"])
lbl <- paste("R^2 == ", 0.963)
monoplot <- ggplot(data = monopred, aes(x=truepropscombined[,"Monocyte"], y=Cell_estimates300[,"Mono"])) + 
  geom_point(color="orange", size=2) + ggtitle("Monoctye Cells") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold"),
        axis.title.x = element_text(size=14, face="bold"), 
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12)) +
  xlab("True Proportions") + ylab("Predicted Proportions") + 
  geom_abline(size=1) + 
  annotate("text",x=10, y=20,label=lbl, parse=TRUE) + 
  annotate("text",x=10, y=18,label="RMSE = 1.55")



#plot for Granulocyte
granpred <- data.frame(Cell_estimates300[,"Gran"], truepropscombined[,"Granulocyte"])
lbl <- paste("R^2 == ", 0.995)
granplot <- ggplot(data = monopred, aes(x=truepropscombined[,"Granulocyte"], y=Cell_estimates300[,"Gran"])) + 
  geom_point(color="deeppink", size=2) + ggtitle("Granulocyte Cells") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold"),
        axis.title.x = element_text(size=14, face="bold"), 
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12)) +
  xlab("True Proportions") + ylab("Predicted Proportions") + 
  geom_abline(size=1) + 
  annotate("text",x=25, y=62,label=lbl, parse=TRUE) + 
  annotate("text",x=25, y=56,label="RMSE = 2.28")

#put all 6 plots on one plot

grid.arrange(cd4plot, cd8plot, bcellplot, nkplot, monoplot, granplot, nrow=2)

################################################################################
#random library

#plot for CD4
cd4predrand <- data.frame(Cell_estimatesrand[,"CD4T"], truepropscombined[,"CD4T"])
lbl <- paste("R^2 == ", 0.825)
cd4plotrand <- ggplot(data = cd4predrand, aes(x=truepropscombined[,"CD4T"], y=Cell_estimatesrand[,"CD4T"])) + 
  geom_point(color="blue", size=2) + ggtitle("CD4") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                                            axis.title.x = element_blank(), 
                                                            axis.title.y = element_blank()) +
  geom_abline(size=1) + 
  annotate("text",x=11, y=25,label=lbl, parse=TRUE) + 
  annotate("text",x=13, y=22,label="RMSE=3.214")



#plot for CD8

cd8predrand <- data.frame(Cell_estimatesrand[,"CD8T"], truepropscombined[,"CD8T"])
lbl <- paste("R^2 == ", 0.965)
cd8plotrand <- ggplot(data = cd8predrand, aes(x=truepropscombined[,"CD8T"], y=Cell_estimatesrand[,"CD8T"])) + 
  geom_point(color="red", size=2) + ggtitle("CD8") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold"),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  geom_abline(size=1) + 
  annotate("text",x=12, y=32,label=lbl, parse=TRUE) + 
  annotate("text",x=14, y=28,label="RMSE=3.056")



#plot for B cell
bcellpredrand <- data.frame(Cell_estimatesrand[,"Bcell"], truepropscombined[,"Bcell"])
lbl <- paste("R^2 == ", 0.971)
bcellplotrand <- ggplot(data = bcellpredrand, aes(x=truepropscombined[,"Bcell"], y=Cell_estimatesrand[,"Bcell"])) + 
  geom_point(color="dark green", size=2) + ggtitle("B") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold"),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  geom_abline(size=1) + 
  annotate("text",x=10, y=27,label=lbl, parse=TRUE) + 
  annotate("text",x=11, y=23,label="RMSE=1.608")



#plot for NK
nkpredrand <- data.frame(Cell_estimatesrand[,"NK"], truepropscombined[,"NK"])
lbl <- paste("R^2 == ", 0.723)
nkplotrand <- ggplot(data = nkpredrand, aes(x=truepropscombined[,"NK"], y=Cell_estimatesrand[,"NK"])) + 
  geom_point(color="purple", size=2) + ggtitle("NK") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold"),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  geom_abline(size=1) + 
  annotate("text",x=7, y=25,label=lbl, parse=TRUE) + 
  annotate("text",x=8, y=21,label="RMSE=4.652")



#plot for Monocytes
monopredrand <- data.frame(Cell_estimatesrand[,"Mono"], truepropscombined[,"Monocyte"])
lbl <- paste("R^2 == ", 0.960)
monoplotrand <- ggplot(data = monopredrand, aes(x=truepropscombined[,"Monocyte"], y=Cell_estimatesrand[,"Mono"])) + 
  geom_point(color="orange", size=2) + ggtitle("Monocyte") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold"),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  geom_abline(size=1) + 
  annotate("text",x=10, y=24,label=lbl, parse=TRUE) + 
  annotate("text",x=11, y=21,label="RMSE=1.686")



#plot for Granulocyte
granpredrand <- data.frame(Cell_estimatesrand[,"Gran"], truepropscombined[,"Granulocyte"])
lbl <- paste("R^2 == ", 0.992)
granplotrand <- ggplot(data = monopredrand, aes(x=truepropscombined[,"Granulocyte"], y=Cell_estimatesrand[,"Gran"])) + 
  geom_point(color="deeppink", size=2) + ggtitle("Granulocyte") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold"),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  geom_abline(size=1) + 
  annotate("text",x=26, y=74,label=lbl, parse=TRUE) + 
  annotate("text",x=31, y=63,label="RMSE= 2.977")

#put all 6 plots on one plot

grid.arrange(cd4plotrand, cd8plotrand, bcellplotrand, 
             nkplotrand, monoplotrand, granplotrand, nrow=2,
             bottom = "True Proportions",
             left = "Predicted Proportions")

###########################################################################
#scale everything to sum to 1 and recalculate RMSE
###########################################################################

scale <- function(x){
  x/sum(x)
}

#center true proportions
truepropscombined_cent <- apply(truepropscombined, 1, scale)
truepropscombined_cent <- t(truepropscombined_cent)

#DSC method
Cell_estimates72_cent <- apply(Cell_estimates72, 1, scale) 
Cell_estimates72_cent <- t(Cell_estimates72_cent)

Cell_estimates120_cent <- apply(Cell_estimates120, 1, scale) 
Cell_estimates120_cent <- t(Cell_estimates120_cent)

Cell_estimates180_cent <- apply(Cell_estimates180, 1, scale) 
Cell_estimates180_cent <- t(Cell_estimates180_cent)

Cell_estimates240_cent <- apply(Cell_estimates240, 1, scale) 
Cell_estimates240_cent <- t(Cell_estimates240_cent)

Cell_estimates300_cent <- apply(Cell_estimates300, 1, scale) 
Cell_estimates300_cent <- t(Cell_estimates300_cent)

Cell_estimates360_cent <- apply(Cell_estimates360, 1, scale) 
Cell_estimates360_cent <- t(Cell_estimates360_cent)

Cell_estimates540_cent <- apply(Cell_estimates540, 1, scale) 
Cell_estimates540_cent <- t(Cell_estimates540_cent)


#RMSE for 72 library
cd4rmse72_cent <- rmse(Cell_estimates72_cent[,"CD4T"], truepropscombined_cent[,"CD4T"])
cd8rmse72_cent <- rmse(Cell_estimates72_cent[,"CD8T"], truepropscombined_cent[,"CD8T"])
nkrmse72_cent <- rmse(Cell_estimates72_cent[,"NK"], truepropscombined_cent[,"NK"])
brmse72_cent <- rmse(Cell_estimates72_cent[,"Bcell"], truepropscombined_cent[,"Bcell"])
monormse72_cent <- rmse(Cell_estimates72_cent[,"Mono"], truepropscombined_cent[,"Monocyte"])
granrmse72_cent <- rmse(Cell_estimates72_cent[,"Gran"], truepropscombined_cent[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse72_cent <- sum(c(cd4rmse72_cent, cd8rmse72_cent, nkrmse72_cent, 
                       brmse72_cent, monormse72_cent, granrmse72_cent))/6

#RMSE for 120 library
cd4rmse120_cent <- rmse(Cell_estimates120_cent[,"CD4T"], truepropscombined_cent[,"CD4T"])
cd8rmse120_cent <- rmse(Cell_estimates120_cent[,"CD8T"], truepropscombined_cent[,"CD8T"])
nkrmse120_cent <- rmse(Cell_estimates120_cent[,"NK"], truepropscombined_cent[,"NK"])
brmse120_cent <- rmse(Cell_estimates120_cent[,"Bcell"], truepropscombined_cent[,"Bcell"])
monormse120_cent <- rmse(Cell_estimates120_cent[,"Mono"], truepropscombined_cent[,"Monocyte"])
granrmse120_cent <- rmse(Cell_estimates120_cent[,"Gran"], truepropscombined_cent[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse120_cent <- sum(c(cd4rmse120_cent, cd8rmse120_cent, nkrmse120_cent, 
                        brmse120_cent, monormse120_cent, granrmse120_cent))/6

#RMSE for 180 library
cd4rmse180_cent <- rmse(Cell_estimates180_cent[,"CD4T"], truepropscombined_cent[,"CD4T"])
cd8rmse180_cent <- rmse(Cell_estimates180_cent[,"CD8T"], truepropscombined_cent[,"CD8T"])
nkrmse180_cent <- rmse(Cell_estimates180_cent[,"NK"], truepropscombined_cent[,"NK"])
brmse180_cent <- rmse(Cell_estimates180_cent[,"Bcell"], truepropscombined_cent[,"Bcell"])
monormse180_cent <- rmse(Cell_estimates180_cent[,"Mono"], truepropscombined_cent[,"Monocyte"])
granrmse180_cent <- rmse(Cell_estimates180_cent[,"Gran"], truepropscombined_cent[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse180_cent <- sum(c(cd4rmse180_cent, cd8rmse180_cent, nkrmse180_cent, 
                        brmse180_cent, monormse180_cent, granrmse180_cent))/6

#RMSE for 240 library
cd4rmse240_cent <- rmse(Cell_estimates240_cent[,"CD4T"], truepropscombined_cent[,"CD4T"])
cd8rmse240_cent <- rmse(Cell_estimates240_cent[,"CD8T"], truepropscombined_cent[,"CD8T"])
nkrmse240_cent <- rmse(Cell_estimates240_cent[,"NK"], truepropscombined_cent[,"NK"])
brmse240_cent <- rmse(Cell_estimates240_cent[,"Bcell"], truepropscombined_cent[,"Bcell"])
monormse240_cent <- rmse(Cell_estimates240_cent[,"Mono"], truepropscombined_cent[,"Monocyte"])
granrmse240_cent <- rmse(Cell_estimates240_cent[,"Gran"], truepropscombined_cent[,6])
#get the average RMSE across all cell types
averagermse240_cent <- sum(c(cd4rmse240_cent, cd8rmse240_cent, nkrmse240_cent, 
                        brmse240_cent, monormse240_cent, granrmse240_cent))/6

#RMSE for 300 library
cd4rmse300_cent <- rmse(Cell_estimates300_cent[,"CD4T"], truepropscombined_cent[,"CD4T"])
cd8rmse300_cent <- rmse(Cell_estimates300_cent[,"CD8T"], truepropscombined_cent[,"CD8T"])
nkrmse300_cent <- rmse(Cell_estimates300_cent[,"NK"], truepropscombined_cent[,"NK"])
brmse300_cent <- rmse(Cell_estimates300_cent[,"Bcell"], truepropscombined_cent[,"Bcell"])
monormse300_cent <- rmse(Cell_estimates300_cent[,"Mono"], truepropscombined_cent[,"Monocyte"])
granrmse300_cent <- rmse(Cell_estimates300_cent[,"Gran"], truepropscombined_cent[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse300_cent <- sum(c(cd4rmse300_cent, cd8rmse300_cent, nkrmse300_cent, 
                        brmse300_cent, monormse300_cent, granrmse300_cent))/6

#RMSE for 360 library
cd4rmse360_cent <- rmse(Cell_estimates360_cent[,"CD4T"], truepropscombined_cent[,"CD4T"])
cd8rmse360_cent <- rmse(Cell_estimates360_cent[,"CD8T"], truepropscombined_cent[,"CD8T"])
nkrmse360_cent <- rmse(Cell_estimates360_cent[,"NK"], truepropscombined_cent[,"NK"])
brmse360_cent <- rmse(Cell_estimates360_cent[,"Bcell"], truepropscombined_cent[,"Bcell"])
monormse360_cent <- rmse(Cell_estimates360_cent[,"Mono"], truepropscombined_cent[,"Monocyte"])
granrmse360_cent <- rmse(Cell_estimates360_cent[,"Gran"], truepropscombined_cent[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse360_cent <- sum(c(cd4rmse360_cent, cd8rmse360_cent, nkrmse360_cent, 
                        brmse360_cent, monormse360_cent, granrmse360_cent))/6

#RMSE for 540 library
cd4rmse540_cent <- rmse(Cell_estimates540_cent[,"CD4T"], truepropscombined_cent[,"CD4T"])
cd8rmse540_cent <- rmse(Cell_estimates540_cent[,"CD8T"], truepropscombined_cent[,"CD8T"])
nkrmse540_cent <- rmse(Cell_estimates540_cent[,"NK"], truepropscombined_cent[,"NK"])
brmse540_cent <- rmse(Cell_estimates540_cent[,"Bcell"], truepropscombined_cent[,"Bcell"])
monormse540_cent <- rmse(Cell_estimates540_cent[,"Mono"], truepropscombined_cent[,"Monocyte"])
granrmse540_cent <- rmse(Cell_estimates540_cent[,"Gran"], truepropscombined_cent[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse540_cent <- sum(c(cd4rmse540_cent, cd8rmse540_cent, nkrmse540_cent, 
                             brmse540_cent, monormse540_cent, granrmse540_cent))/6

centered_rmse <- matrix(NA, nrow = 6, ncol = 7)
rownames(centered_rmse) <- c("CD4", "CD8", "NK","B","Mono","Gran")
colnames(centered_rmse) <- c("72","120","180","240","300","360","540")
centered_rmse[1,] <- c(cd4rmse72_cent, cd4rmse120_cent, cd4rmse180_cent,cd4rmse240_cent,
                       cd4rmse300_cent,cd4rmse360_cent,cd4rmse540_cent)
centered_rmse[2,] <- c(cd8rmse72_cent,cd8rmse120_cent,cd8rmse180_cent,cd8rmse240_cent,
                       cd8rmse300_cent,cd8rmse360_cent,cd8rmse540_cent)
centered_rmse[3,] <- c(nkrmse72_cent, nkrmse120_cent,nkrmse180_cent,nkrmse240_cent,
                       nkrmse300_cent,nkrmse360_cent,nkrmse540_cent)
centered_rmse[4,] <- c(brmse72_cent, brmse120_cent,brmse180_cent,brmse240_cent,
                       brmse300_cent,brmse360_cent,brmse540_cent)
centered_rmse[5,] <- c(monormse72_cent, monormse120_cent,monormse180_cent,monormse240_cent,
                       monormse300_cent,monormse360_cent,monormse540_cent)
centered_rmse[6,] <- c(granrmse72_cent,granrmse120_cent,granrmse180_cent,granrmse240_cent,
                       granrmse300_cent,granrmse360_cent,granrmse540_cent)

#legacy method
Cell_estimates72_t_cent <- apply(Cell_estimates72_t, 1, scale) 
Cell_estimates72_t_cent <- t(Cell_estimates72_t_cent)

Cell_estimates120_t_cent <- apply(Cell_estimates120_t, 1, scale) 
Cell_estimates120_t_cent <- t(Cell_estimates120_t_cent)

Cell_estimates180_t_cent <- apply(Cell_estimates180_t, 1, scale) 
Cell_estimates180_t_cent <- t(Cell_estimates180_t_cent)

Cell_estimates240_t_cent <- apply(Cell_estimates240_t, 1, scale) 
Cell_estimates240_t_cent <- t(Cell_estimates240_t_cent)

Cell_estimates300_t_cent <- apply(Cell_estimates300_t, 1, scale) 
Cell_estimates300_t_cent <- t(Cell_estimates300_t_cent)

Cell_estimates360_t_cent <- apply(Cell_estimates360_t, 1, scale) 
Cell_estimates360_t_cent <- t(Cell_estimates360_t_cent)

Cell_estimates540_t_cent <- apply(Cell_estimates540_t, 1, scale) 
Cell_estimates540_t_cent <- t(Cell_estimates540_t_cent)


#now get rmse for 72 library
cd4rmse72_t_cent <- rmse(Cell_estimates72_t_cent[,"CD4T"], truepropscombined_cent[,"CD4T"])
cd8rmse72_t_cent <- rmse(Cell_estimates72_t_cent[,"CD8T"], truepropscombined_cent[,"CD8T"])
nkrmse72_t_cent <- rmse(Cell_estimates72_t_cent[,"NK"], truepropscombined_cent[,"NK"])
brmse72_t_cent <- rmse(Cell_estimates72_t_cent[,"Bcell"], truepropscombined_cent[,"Bcell"])
monormse72_t_cent <- rmse(Cell_estimates72_t_cent[,"Mono"], truepropscombined_cent[,"Monocyte"])
granrmse72_t_cent <- rmse(Cell_estimates72_t_cent[,"Gran"], truepropscombined_cent[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse72_t_cent <- sum(c(cd4rmse72_t_cent, cd8rmse72_t_cent, nkrmse72_t_cent,
                         brmse72_t_cent, monormse72_t_cent, granrmse72_t_cent))/6

#now get rmse for 120 library
cd4rmse120_t_cent <- rmse(Cell_estimates120_t_cent[,"CD4T"], truepropscombined_cent[,"CD4T"])
cd8rmse120_t_cent <- rmse(Cell_estimates120_t_cent[,"CD8T"], truepropscombined_cent[,"CD8T"])
nkrmse120_t_cent <- rmse(Cell_estimates120_t_cent[,"NK"], truepropscombined_cent[,"NK"])
brmse120_t_cent <- rmse(Cell_estimates120_t_cent[,"Bcell"], truepropscombined_cent[,"Bcell"])
monormse120_t_cent <- rmse(Cell_estimates120_t_cent[,"Mono"], truepropscombined_cent[,"Monocyte"])
granrmse120_t_cent <- rmse(Cell_estimates120_t_cent[,"Gran"], truepropscombined_cent[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse120_t_cent <- sum(c(cd4rmse120_t_cent, cd8rmse120_t_cent, nkrmse120_t_cent, 
                          brmse120_t_cent, monormse120_t_cent, granrmse120_t_cent))/6

#now get rmse for 180 library
cd4rmse180_t_cent <- rmse(Cell_estimates180_t_cent[,"CD4T"], truepropscombined_cent[,"CD4T"])
cd8rmse180_t_cent <- rmse(Cell_estimates180_t_cent[,"CD8T"], truepropscombined_cent[,"CD8T"])
nkrmse180_t_cent <- rmse(Cell_estimates180_t_cent[,"NK"], truepropscombined_cent[,"NK"])
brmse180_t_cent <- rmse(Cell_estimates180_t_cent[,"Bcell"], truepropscombined_cent[,"Bcell"])
monormse180_t_cent <- rmse(Cell_estimates180_t_cent[,"Mono"], truepropscombined_cent[,"Monocyte"])
granrmse180_t_cent <- rmse(Cell_estimates180_t_cent[,"Gran"], truepropscombined_cent[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse180_t_cent <- sum(c(cd4rmse180_t_cent, cd8rmse180_t_cent, nkrmse180_t_cent, 
                          brmse180_t_cent, monormse180_t_cent, granrmse180_t_cent))/6

#now get rmse for 240 library
cd4rmse240_t_cent <- rmse(Cell_estimates240_t_cent[,"CD4T"], truepropscombined_cent[,"CD4T"])
cd8rmse240_t_cent <- rmse(Cell_estimates240_t_cent[,"CD8T"], truepropscombined_cent[,"CD8T"])
nkrmse240_t_cent <- rmse(Cell_estimates240_t_cent[,"NK"], truepropscombined_cent[,"NK"])
brmse240_t_cent <- rmse(Cell_estimates240_t_cent[,"Bcell"], truepropscombined_cent[,"Bcell"])
monormse240_t_cent <- rmse(Cell_estimates240_t_cent[,"Mono"], truepropscombined_cent[,"Monocyte"])
granrmse240_t_cent <- rmse(Cell_estimates240_t_cent[,"Gran"], truepropscombined_cent[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse240_t_cent <- sum(c(cd4rmse240_t_cent, cd8rmse240_t_cent, nkrmse240_t_cent,
                          brmse240_t_cent, monormse240_t_cent, granrmse240_t_cent))/6

#now get rmse for 300 library
cd4rmse300_t_cent <- rmse(Cell_estimates300_t_cent[,"CD4T"], truepropscombined_cent[,"CD4T"])
cd8rmse300_t_cent <- rmse(Cell_estimates300_t_cent[,"CD8T"], truepropscombined_cent[,"CD8T"])
nkrmse300_t_cent <- rmse(Cell_estimates300_t_cent[,"NK"], truepropscombined_cent[,"NK"])
brmse300_t_cent <- rmse(Cell_estimates300_t_cent[,"Bcell"], truepropscombined_cent[,"Bcell"])
monormse300_t_cent <- rmse(Cell_estimates300_t_cent[,"Mono"], truepropscombined_cent[,"Monocyte"])
granrmse300_t_cent <- rmse(Cell_estimates300_t_cent[,"Gran"], truepropscombined_cent[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse300_t_cent <- sum(c(cd4rmse300_t_cent, cd8rmse300_t_cent, nkrmse300_t_cent,
                          brmse300_t_cent, monormse300_t_cent, granrmse300_t_cent))/6

#now get rmse for 360 library
cd4rmse360_t_cent <- rmse(Cell_estimates360_t_cent[,"CD4T"], truepropscombined_cent[,"CD4T"])
cd8rmse360_t_cent <- rmse(Cell_estimates360_t_cent[,"CD8T"], truepropscombined_cent[,"CD8T"])
nkrmse360_t_cent <- rmse(Cell_estimates360_t_cent[,"NK"], truepropscombined_cent[,"NK"])
brmse360_t_cent <- rmse(Cell_estimates360_t_cent[,"Bcell"], truepropscombined_cent[,"Bcell"])
monormse360_t_cent <- rmse(Cell_estimates360_t_cent[,"Mono"], truepropscombined_cent[,"Monocyte"])
granrmse360_t_cent <- rmse(Cell_estimates360_t_cent[,"Gran"], truepropscombined_cent[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse360_t_cent <- sum(c(cd4rmse360_t_cent, cd8rmse360_t_cent, nkrmse360_t_cent,
                          brmse360_t_cent, monormse360_t_cent, granrmse360_t_cent))/6

#now get rmse for 540 library
cd4rmse540_t_cent <- rmse(Cell_estimates540_t_cent[,"CD4T"], truepropscombined_cent[,"CD4T"])
cd8rmse540_t_cent <- rmse(Cell_estimates540_t_cent[,"CD8T"], truepropscombined_cent[,"CD8T"])
nkrmse540_t_cent <- rmse(Cell_estimates540_t_cent[,"NK"], truepropscombined_cent[,"NK"])
brmse540_t_cent <- rmse(Cell_estimates540_t_cent[,"Bcell"], truepropscombined_cent[,"Bcell"])
monormse540_t_cent <- rmse(Cell_estimates540_t_cent[,"Mono"], truepropscombined_cent[,"Monocyte"])
granrmse540_t_cent <- rmse(Cell_estimates540_t_cent[,"Gran"], truepropscombined_cent[,"Granulocyte"])
#get the average RMSE across all cell types
averagermse540_t_cent <- sum(c(cd4rmse540_t_cent, cd8rmse540_t_cent, nkrmse540_t_cent,
                          brmse540_t_cent, monormse540_t_cent, granrmse540_t_cent))/6


centered_rmse_t <- matrix(NA, nrow = 6, ncol = 7)
rownames(centered_rmse_t) <- c("CD4", "CD8", "NK","B","Mono","Gran")
colnames(centered_rmse_t) <- c("72","120","180","240","300","360","540")
centered_rmse_t[1,] <- c(cd4rmse72_t_cent, cd4rmse120_t_cent, cd4rmse180_t_cent,cd4rmse240_t_cent,
                       cd4rmse300_t_cent,cd4rmse360_t_cent,cd4rmse540_t_cent)
centered_rmse_t[2,] <- c(cd8rmse72_t_cent,cd8rmse120_t_cent,cd8rmse180_t_cent,cd8rmse240_t_cent,
                       cd8rmse300_t_cent,cd8rmse360_t_cent,cd8rmse540_t_cent)
centered_rmse_t[3,] <- c(nkrmse72_t_cent, nkrmse120_t_cent,nkrmse180_t_cent,nkrmse240_t_cent,
                       nkrmse300_t_cent,nkrmse360_t_cent,nkrmse540_t_cent)
centered_rmse_t[4,] <- c(brmse72_t_cent, brmse120_t_cent,brmse180_t_cent,brmse240_t_cent,
                       brmse300_t_cent,brmse360_t_cent,brmse540_t_cent)
centered_rmse_t[5,] <- c(monormse72_t_cent, monormse120_t_cent,monormse180_t_cent,monormse240_t_cent,
                       monormse300_t_cent,monormse360_t_cent,monormse540_t_cent)
centered_rmse_t[6,] <- c(granrmse72_t_cent,granrmse120_t_cent,granrmse180_t_cent,granrmse240_t_cent,
                       granrmse300_t_cent,granrmse360_t_cent,granrmse540_t_cent)

####now square and average over library sizes

centered_rmsesq <- centered_rmse^2
apply(centered_rmsesq, 1, FUN = mean)

centered_rmse_tsq <- centered_rmse_t^2
apply(centered_rmse_tsq, 1, FUN = mean)




