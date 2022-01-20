#load packages
library(quadprog)

#load functions and libraries
load("Reference Data for IDOL.RData")  # reference data
source("UnsupervisedOptimalDMRFinderV3.R")
source("IDOL revised functions 02_19_2017.R")
source("multiRegMultistat.R")
load("Legacy240.RData")
load("Optimal set of 240 CpGs.RData")


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
load("HannumFullAnalysisWorkspace.RData")

###########################
#Using modified DSC library
###########################

# cell predictions
interDSCHan = intersect(library240, rownames(meth.data))

HannumDSC = meth.data[interDSCHan,]
HanmeanmethDSC = meanmeth[interDSCHan,]

cellTypes = colnames(HanmeanmethDSC)
Lwbc = diag(6)
colnames(Lwbc) = rownames(Lwbc) = cellTypes
CellPredsDSCHan240 = data.frame(projectWBCnew(HannumDSC, HanmeanmethDSC, Lwbc)*100)
save(CellPredsDSCHan240, file= "Cell Preds DSC Hannum 240.RData")

# compute the multivariate R2
meth.data.sub = meth.data
XDSCHan = model.matrix(~ CD8T + CD4T + NK + Bcell + Gran, data = CellPredsDSCHan240)

Llist = list()
Llist[[1]] = c(1,0,0,0,0,0)
MultiDSCHan240 = multiRegMultistat(meth.data.sub, XDSCHan, Llist)

save(MultiDSCHan240, file = "Hannum modifed DSC Results 240.RData")

###########################
#Using Legacy library
###########################

# cell predictions

interlegHan = intersect(tstats240, rownames(meth.data))

HannumLeg = meth.data[interlegHan,]
Hanmeanmethleg = meanmeth[interlegHan,]

cellTypes = colnames(Hanmeanmethleg)
Lwbc = diag(6)
colnames(Lwbc) = rownames(Lwbc) = cellTypes
CellPredslegHan240 = data.frame(projectWBCnew(HannumLeg, Hanmeanmethleg, Lwbc)*100)
save(CellPredslegHan240, file= "Cell Preds Legacy Hannum 240.RData")

# compute the multivariate R2
meth.data.sub = meth.data
Xleghan = model.matrix(~ CD8T + CD4T + NK + Bcell + Gran, data = CellPredslegHan240)

Llist = list()
Llist[[1]] = c(1,0,0,0,0,0)
MultilegHan240 = multiRegMultistat(meth.data.sub, Xleghan, Llist)

save(MultilegHan240, file = "Hannum Legacy Results 240.RData")