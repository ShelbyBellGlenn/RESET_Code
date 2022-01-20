#load packages
library(quadprog)

#load functions and libraries
load("Reference Data for IDOL.RData")  # reference data
source("UnsupervisedOptimalDMRFinderV3.R")
source("IDOL revised functions 02_19_2017.R")
source("multiRegMultistat.R")
load("Legacy360.RData")
load("Optimal set of 360 CpGs.RData")


#get mean methylation matrix, calculate from the reference betas matrix
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


#########################################################################
#Feinberg analysis with library size of 120
#########################################################################
load("ObjectsForFeinberg.RData")

###########################
#Using modified DSC library
###########################

#subset feinberg data to our library
DSCOptimFein <- dat[library360,]

#subset mean methylation matrix to our library
DSCmeanmethFein <- meanmeth[library360,]

#get cell proportion estimates
cellTypes = colnames(DSCmeanmethFein)
Lwbc = diag(6)
colnames(Lwbc) = rownames(Lwbc) = cellTypes
CellPredsDSCFein360 = data.frame(projectWBCnew(DSCOptimFein, DSCmeanmethFein , Lwbc)*100)
save(CellPredsDSCFein360, file= "Cell Preds DSC Feinberg 360.RData")

#get design matrix for regression
XDSCFein = model.matrix(~ CD8T + CD4T + NK + Bcell + Gran, data = CellPredsDSCFein360)

#compute the multivariate R2
Llist = list()
Llist[[1]] = c(1,0,0,0,0,0)
MultiDSCFein360 = multiRegMultistat(dat, XDSCFein, Llist)

save(MultiDSCFein360, file = "Feinberg Modified DSC Results 360.RData")

###########################
#Using Legacy library
###########################

#subset feinberg data to our library
legacyOptimFein <- dat[tstats360,]

#subset mean methylation matrix to our library
legmeanmethFein <- meanmeth[tstats360,]

#get cell proportion estimates
cellTypes = colnames(legmeanmethFein)
Lwbc = diag(6)
colnames(Lwbc) = rownames(Lwbc) = cellTypes
CellPredslegFein360 = data.frame(projectWBCnew(legacyOptimFein, legmeanmethFein, Lwbc)*100)
save(CellPredslegFein360, file= "Cell Preds Legacy Feinberg 360.RData")

#get design matrix for regression
XLegFein = model.matrix(~ CD8T + CD4T + NK + Bcell + Gran, data = CellPredslegFein360)

# compute the multivariate R2
Llist = list()
Llist[[1]] = c(1,0,0,0,0,0)
MultiLegFein360 = multiRegMultistat(dat, XLegFein, Llist)

save(MultiLegFein360, file = "Feinberg Legacy Results 360.RData")
