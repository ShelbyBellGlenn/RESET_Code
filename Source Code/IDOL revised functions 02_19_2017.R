#############################################################################################
#############################################################################################
#  Author: Devin C. Koestler
#  Date: February 19, 2017
#  Title: Updated functions for IDOL (Koestler et al., 2016)
#############################################################################################
#############################################################################################

#############################################################################################
# FUNCTION:  CandidateDMRFinder.v2 
#    This function identifies candidate/putative differentially methylated loci (DML)
#    based on the procedure described in Koestler et al., (2016).  Breifly, a series
#    of two-sample t-tests are fit to the J CpGs contained in the referenceBetas 
#    object and used to compare the mean methylation beta-values between each of the
#    K cell type against the mean methylation beta-values computed across the remaining 
#    K - 1 cell types.  Putative DMLs are identified by first rank ordering CpGs by their 
#    t-statistics, then taking the top M DMLs with the smallest and largest t-statistics 
#    for each of the K comparisons.  
#
# REQUIRES:	genefilter    
#
# ARGUMENTS:
#   cellTypes:       A vector of length K that contains the cell type names.  For example,
#					 c("CD4T", "CD8T", "NK", "Bcell", "Mono", "Gran").
#                    
#	referenceBetas:  A J x N matrix of cell-specific methylation beta-values; J represents 
#                    the number of CpGs (i.e., ~ 450,000 for the Illumina HumanMethylation450
#                    array) and N represents the number of samples for which cell-specific
#                    methylation signatures are available.
#
#
#	referenceCovars: A N x P data.frame of meta data aross the N samples.  The rows of this
#                    object MUST be in the same order as the columns of referenceBetas.
#                    Further, there must be a column called "CellType" (case sensitive), 
#                    that indicates the cell-type identity for each of the N samples.
#                    The nomenclature used to indicate cell identity across the N samples
#                    should follow the nomenclature used for cellTypes (see above).
#
#	M:				 The number of candidate DMLs with the smallest and largest t-statistic
#                    to return for each comparison.  Defaults to M = 150 as in Koestler et al.,
#                    (2016)
#
#	equal.variance:  Should a t-test assuming equal variances be fit.  Defaults to FALSE, 
#                    an unequal variance t-test.
#
# RETURNS:   A list containing two objects: (1) candidateSet - a vector containing the names 
#            of the R candidate DMLs identified from the analysis and (2) coefEsts - A R x K
#            matrix of the within-cell type mean methylation beta values across the R 
#            identified candidate DMLs.
#    
#############################################################################################

CandidateDMRFinder.v2 =function(cellTypes, referenceBetas, referenceCovars, M = 150, equal.variance = F){

	require(genefilter)

	p = referenceBetas
	pd = referenceCovars
	K = length(cellTypes)

	if(sum(cellTypes %in% pd$CellType)!= K) {
		stop("cell type names in target covariate data are not consistent with the nomenclature used in cellTypes")
	}

	splitit <- function(x) {
        split(seq(along = x), x)
    }

    keep <- which(pd$CellType %in% cellTypes)
    pd <- pd[keep, ]
    p <- p[, keep]

    tIndexes <- splitit(pd$CellType)
    tstatList1 <- lapply(tIndexes, function(i) {
       x1 <- i
       x2 <- c(1:ncol(p))[-x1]
       return(fastT(p, x1, x2, var.equal = equal.variance))
    })

    probeList1 <- lapply(tstatList1, function(x) {
        yUp <- rownames(p)[order(x$z, decreasing = TRUE)]
        yDown <- rownames(p)[order(x$z, decreasing = FALSE)]
        c(yUp[1:M], yDown[1:M])
    })

    candidateSet <- unique(unlist(probeList1))
    p <- p[candidateSet, ]

    coefEsts = matrix(NA, nrow = length(candidateSet), ncol = K)
    rownames(coefEsts) = candidateSet
    colnames(coefEsts) = cellTypes
    for(k in 1:K) {
    	ind = which(pd$CellType %in% cellTypes[k])
    	coefEsts[,k] = apply(p[,ind], 1, mean, na.rm = T)
    }

    tmp = list(candidateSet, coefEsts)
    names(tmp) = c("candidateSet", "coefEsts")
    tmp
}

#############################################################################################
# FUNCTION:  projectWBCnew
#    This function predicts the underlying cellular composition of heterogeneous tissue 
#    types (i.e., WB) using the constrained projection procedure described by Houseman et al., 
#    (2012).  
#
# REQUIRES:	quadprog    
#
# ARGUMENTS:
#      
#	Y:  			 A J x N matrix of methylation beta-values collected from mixed/
#                    heterogeneous biospecimen (i.e., WB).  Target set.
#
#	coefWBC:         A J x K projection matrix;, i.e., within-cell type mean methylation 
#                    matrix across J DMLs and K many cell types
#
#	contrastWBC:	 Contrast for cell composition predictions.  The user needn't modify
#                    this 
#
#	nonnegative:     Should cell predictions be nonnegative.  Defaults to TRUE
#
# RETURNS:   A N x K matrix of cell proportion estimates across the K cell types for each 
#            of the N subjects contained in the Target Set.
#    
#############################################################################################

projectWBCnew = function(Y, coefWBC, contrastWBC=NULL, nonnegative=TRUE){ 
	
  if(is.null(contrastWBC)) Xmat = coefWBC
  else Xmat = coefWBC %*% t(contrastWBC) 

  nCol = dim(Xmat)[2]
  nSubj = dim(Y)[2]

  mixCoef = matrix(0, nSubj, nCol)
  rownames(mixCoef) = colnames(Y)
  colnames(mixCoef) = colnames(Xmat)

  if(nonnegative){
    library(quadprog)

   Amat = cbind(rep(-1,nCol), diag(nCol))
   b0vec = c(-1,rep(0,nCol))

    for(i in 1:nSubj){
      obs = which(!is.na(Y[,i])) 
      Dmat = t(Xmat[obs,])%*%Xmat[obs,]
      mixCoef[i,] = solve.QP(Dmat, t(Xmat[obs,])%*%Y[obs,i], Amat, b0vec, meq = 0)$sol
    }
  }
  else{
    for(i in 1:nSubj){
      obs = which(!is.na(Y[,i])) 
      Dmat = t(Xmat[obs,])%*%Xmat[obs,]
      mixCoef[i,] = solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
    }
  }

  return(mixCoef)
}

#############################################################################################
# FUNCTION:  IDOLoptimize
#    This function identifies optimal DML/DMR libraries for cell mixture deconvolution using 
#    the procedure described in Koestler et al., (2016).  
#
# REQUIRES:	quadprog and doParallel     
#
# ARGUMENTS:
#   candDMRFinderObject:    List object returned from the CandidateDMRFinder.v2 function 
#                    
#	trainingBetas:  		A J x N matrix of beta-values where methylation was profiled 
#                           in a heterogenous tissue type (i.e., WB).  Here, J indicates
#                           the number of CpGs and N, the number of samples
#
#	trainingCovariates: 	A N x P data.frame of meta data aross the N samples.  Contained
#                           within the meta data must be the observed cell fractions (i.e. FACS or
#                           otherwise) for the K cell types indicated in the list object returned 
#                           from the CandidateDMRFinder.v2 function.  The naming of cell types
#                           should be consistent between these two objects as well.
#                           
#	libSize:				Size of the optimized IDOL library.  Defaults to 300.
#
#	maxIt:  				Maximum number of iterations for the IDOL algorithm. Defaults 
#                           to 500.
#
#   numCores:               Number of processing cores to use (see R package DoParallel).  
#                           Defaults to 4.
#
# RETURNS:   A list containing six objects: (1) "IDOL Optimized Library" - a vector containing the names 
#            of the CpGs in the identified IDOL optimized library (2) "IDOL Optimized CoefEsts" - 
#            matrix of the within-cell type mean methylation beta values across the CpGs in the
#            optimal library (3) "RMSE" average root-mean squared error calculated across each 
#            iteration of IDOL (4) "R2" average R2 (coefficent of determinatino) calculated across each 
#            iteration of IDOL (5) "Number of iterations" how many iterations of IDOL were used
#            (6) "Library Size" libsize above.
# 
#############################################################################################

IDOLoptimize = function(candDMRFinderObject, trainingBetas, trainingCovariates, 
	                    libSize = 300, maxIt = 500, numCores = 4) {
    
	# Load necessary packages
	require(quadprog)
	require(doParallel)
	registerDoParallel(cores=numCores)

	# Define relevant functions
	expit = function(w) exp(w)/(1 + exp(w))
	logit = function(w) log(w) - log(1-w)

	R2compute = function(obs, pred) {
		r2 = NULL
		for(i in 1:dim(obs)[2]) {
			y = obs[,i]
			x = pred[,i]
			r2[i] = summary(lm(y~x))[[8]]
		}
		sum(r2)/dim(obs)[2]
	}

	# Polar coordinates function 
	polar = function(x, y, scale = 1) {
		r = sqrt(x^2 + y^2)
		theta = atan2(y, x)
		r*cos(theta - (scale*pi/4))
	}

	# Define relevant parameters
	trainingProbes1 = candDMRFinderObject$candidateSet
	coefEsts = candDMRFinderObject$coefEsts
	P = length(trainingProbes1)
	ProbVector = rep(1/P, P)
	V = libSize  # number of CpGs to select
	B = maxIt # number of interations

	R2.null = 0
	R2Vals = NULL
	RMSE.null = 10000
	RMSEVals = NULL

	# do some cross-checks before beginning algorithm
	cellTypes = colnames(coefEsts)
	K = length(cellTypes)
	if(sum(cellTypes %in% colnames(trainingCovariates)) != K) {
		stop("cell type names in target covariate data are not consistent with cell types to be deconvoluted")
	}

	# Perform the idol maximization algorithm
	for(i in 1:B) {
		Probes = sample(1:P, V, prob = ProbVector)
		CpGNames = trainingProbes1[Probes]

		Beta = coefEsts[CpGNames, ]
		Lwbc = diag(ncol(Beta))

		ctpred = data.frame(projectWBCnew(trainingBetas[CpGNames,], Beta))
		omega.tilde = 100*ctpred
		omega.obs = trainingCovariates[,cellTypes]

		# compute RMSE based on all probes
		# Whole blood
		diff = as.matrix(omega.tilde - omega.obs)
		RMSE = sqrt(sum(diff^2)/K)
		
		# compute R2 based on all probes
		R2 = R2compute(omega.obs, omega.tilde)

		# compute the leave-one-out R-squared values for each of the probes
		Perform.q = foreach(j = 1:length(CpGNames)) %dopar% {

			Beta.q = Beta[CpGNames[-j],]
		
			ctpred.q = data.frame(projectWBCnew(trainingBetas[CpGNames[-j],], Beta.q, Lwbc))
			omega.tilde.q = 100*ctpred.q

			# compute RMSE based on all probes
			# Whole blood
			diff.q = as.matrix(omega.tilde.q - omega.obs)
			RMSE.q = sqrt(sum(diff.q^2)/K)
		
			# compute R2 based on all probes
			R2.q = R2compute(omega.obs, omega.tilde.q)
	    	cbind(R2.q, RMSE.q)
		}
		
		R2.q = unlist(Perform.q)[seq(from = 1, to = (2*V-1), by = 2)]
		RMSE.q = unlist(Perform.q)[seq(from = 2, to = 2*V, by = 2)]

		# compute relevant values to modify the probabilities that a probe is selected in subsequent iterations
		rmse.dq = (RMSE - RMSE.q)*(-1)
		norm.rmse = (rmse.dq)/sd(rmse.dq)
		
		r2.dq = (R2 - R2.q)
		norm.r2 = (r2.dq)/sd(r2.dq)
		
		p1 = polar(norm.rmse, norm.r2)
		
		for (j in 1:length(Probes)) {
			p0 = ProbVector[[Probes[j]]]
			ProbVector[[Probes[j]]] = expit(p1[j])*p0 + p0/2
		}

		# rescale the ProbMat
		ProbVector = ProbVector/sum(ProbVector)

		# save the optimal library for later use
		if(RMSE <= RMSE.null & R2 >= R2.null) {
			RMSE.null = RMSE
			R2.null = R2
			print(paste("Iteration: ", i, " RMSE=", round(RMSE, 3), "; R2=", round(R2, 3), sep = ""))
			
			IDOL.optim.DMRs = CpGNames
			IDOL.optim.coefEsts = coefEsts[CpGNames,]

			save(IDOL.optim.DMRs, IDOL.optim.coefEsts, file = paste("IDOL optimized DMR library_", V, ".RData", sep = ""))
		}
		RMSEVals[i] = RMSE
		R2Vals[i] = R2
	}
	IDOLObjects = list(IDOL.optim.DMRs, IDOL.optim.coefEsts, RMSEVals, R2Vals, B, V)
	names(IDOLObjects) = c("IDOL Optimized Library", "IDOL Optimized CoefEsts", 
		                   "RMSE", "R2", "Number of Iterations", "LibrarySize")

	print(paste("The Average RMSE = ", round(RMSE.null, 3), " and R2 = ", round(R2.null, 3), 
		        " for the IDOL Optimized Library", sep = ""))

	return(IDOLObjects)
}





