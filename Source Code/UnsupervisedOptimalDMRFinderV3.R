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
# RETURNS:   A list containing two objects: (1) candidateDMRs - a vector containing the names 
#            of the R candidate DMLs identified from the analysis and (2) meanMeth - A R x K
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
    names(tmp) = c("candidateDMRs", "meanMeth")
    tmp
}

#########################################################################
# FUNCTION:   DSC.v2
#########################################################################
#########################################################################
# ARGUMENTS:
#
# 	RefBetasSub:    A U x n matrix of methylation beta-values for samples 
#                   where DNAm was profiled in isolated leukocyte subtypes 
#                   U indicates the number of candidate DMRs and n
#                   indicates the number of samples 
#
# 	CellTypes:      A vector of length n whose elements are the cell 
#                   identities of the corresponding to the columns
#                   of RefBetasSub 
#########################################################################
# OUTPUT:        	A vector of DSC values whose first element is the 
#                  overall DSC and whose following elements represent
#                    the DSCs computed between each pair of cell types
#########################################################################

DSC.v2 <- function(RefBetasSub, CellTypes) {
  data = t(RefBetasSub)
	n = nrow(data)
	J = ncol(data)
	batchlevs = levels(factor(CellTypes))
	numPairs = choose(length(batchlevs), 2)
	Pairs = combn(batchlevs, 2)
	PairNames = apply(Pairs, 2, function(w) paste(w[1], " vs ", w[2], sep = ""))

	overallweights = table(factor(CellTypes))/length(CellTypes)
	weights = list()
	for(p in 1:(numPairs+1)) {
		if(p == 1) weights[[p]] = overallweights
		else {
			tmp = rep(0, length(batchlevs))
			tmp[batchlevs %in% Pairs[,(p-1)]] = 0.5
			names(tmp) = batchlevs
			weights[[p]] = tmp
		}
	}
    names(weights) = c("Overall", PairNames)

	avg.w = NULL
	cent.within = matrix(nrow = J, ncol = length(batchlevs))
	rownames(cent.within) = colnames(data)
	colnames(cent.within) = batchlevs
	for(i in 1:length(batchlevs)) {
		dat.w = as.matrix(data[CellTypes == batchlevs[i],])
    if (ncol(dat.w) == 1) dat.w = t(dat.w)
		n.sub = nrow(dat.w)
		cent.w = t(dat.w) %*% rep(1/n.sub, n.sub)
		cent.within[,i] = as.vector(cent.w)

		dist.w = foreach(j=1:nrow(dat.w)) %do% {
			x = as.numeric(dat.w[j,])
			d = x-cent.w
			sqrt(t(d) %*% d)
		}
		avg.w[i] = mean(unlist(dist.w), na.rm = T)
	}
	names(avg.w) = batchlevs

	DSCvals = lapply(weights, function(r) { 
        wt = r
        cells = batchlevs[wt!=0]
        samp.ind = CellTypes %in% cells
        data.ind = data[samp.ind,]
        n.ind = nrow(data.ind)
        cent.global.ind = t(data.ind) %*% rep(1/n.ind, n.ind)
        cent.within.ind = cent.within[,cells]
        dist.b = apply(cent.within.ind, 2, function(v) {
            dist.ind = v - as.vector(cent.global.ind)
            sqrt(t(dist.ind) %*% dist.ind)
            })
		Db = wt[cells] %*% dist.b
		Dw = wt[cells] %*% avg.w[cells]
		DSC = Db/Dw
    })	

    DSCvals = unlist(DSCvals)	
    return(DSCvals)
}

##########################################################################
# FUNCTION:   DMRsubsetFinder
##########################################################################
#
# ARGUMENTS:
#
#   candidateDMRs:  A vector of U candidate DMRs from which DMR subsets should
#                   be randomly selected from, i.e., the output from step 1.
#                
# 	RefBetasSub:    A U x n matrix of methylation beta-values for samples 
#                   where DNAm was profiled in isolated leukocyte subtypes 
#                   U indicates the number of candidate DMRs and n
#                   indicates the number of samples 
#
# 	CellTypes:      A vector of length n whose elements are the cell 
#                   identities of the corresponding to the columns
#                   of RefBetasSub
#
# 	Leuko:       	A vector of length K whose elements are the names of 
#                	of the unique cell types contained in CellTypes, i.e., 
#                	Leuko <- unique(CellTypes)  
#
#   numIter         Number of randomly selected DMR subsets for which 
#                   the DSC is computed on.  Defaults to 1000.
#
# 	numDMRs:     	Vector of length 2 whose elements are the minimum and 
#                   maximum number of DMRs to be selected from the 
#                   candidate DMR search space, i.e., candidateDMRs. Defaults
#                   to a minimum of 50 and a maximum of U-50.
#                      
# OUTPUT:        	A list that is comprised of the following three objects:
#
#	DMRsubsets: 	A list of length numIter whose elements represent vectors
#                   of the CpG names for the numIter randomly selected DMR subsets
# 
#   subsetSizes:    A vector of length numIter whose elements represent 
#                   the size of each randomly selected DMR subset.  Corresponds
#                   to the length of each element of the DMRsubsets list object
#                     
#   DSCvalues:      A list of length numIter whose elements represent the 
#                   computed DSC values based on each of the numIter randomly
#                   selected DMR subsets.
#########################################################################

DMRsubsetFinder <- function(candidateDMRs, RefBetasSub, CellTypes, Leuko,
                   numIter = 1000, numDMRs = c(50, (length(candidateDMRs)-50))) {

	if(!is.matrix(RefBetasSub)) {
    	RefBetasSub = as.matrix(RefBetasSub)
    }  
    
    #########################################################################
    # Randomly select DMR subsets of random sizes from candidateDMRs
    ######################################################################### 
    minSize = numDMRs[1]
    maxSize = numDMRs[2]
    DMRsubsets = list()
    for(j in 1:numIter) {
    	Jstar = round(runif(1, min = minSize, max = maxSize))	   
    	subsetMembers = sample(candidateDMRs, size = Jstar, replace = FALSE) 
    	DMRsubsets[[j]] = subsetMembers
    }
    names(DMRsubsets) = paste("DMRsubset_", 1:numIter, sep = "")
    subsetSize = unlist(lapply(DMRsubsets, length))
    names(subsetSize) = names(DMRsubsets) 

    #########################################################################
    # Compute DSC values for each of the randomly selected DMR subsets
    #########################################################################
    DSCvals = lapply(DMRsubsets, function(b) {
    	Probes = b
    	BetasDMRSubset = RefBetasSub[Probes,]
    	DSC.v2(BetasDMRSubset, CellTypes)
    	})
    names(DSCvals) = names(DMRsubsets)

    #########################################################################
    # Output the results
    #########################################################################
    out <- list(DMRsubsets = DMRsubsets, subsetSizes = subsetSize, 
                DSCvalues = DSCvals)

    return(out) 
}          


####################################################################################
#FUNCTION:   NewDSC.v1
####################################################################################
# ARGUMENTS:
#
# 	RefBetasSub:    A U x n matrix of methylation beta-values for samples 
#                   where DNAm was profiled in isolated leukocyte subtypes 
#                   U indicates the number of candidate DMRs and n
#                   indicates the number of samples 
#
# 	CellTypes:      A vector of length n whose elements are the cell 
#                   identities of the corresponding to the columns
#                   of RefBetasSub 
#########################################################################
# OUTPUT:        	The new overall DSC value computed
#                  
#                    
#########################################################################

NewDSC.v1 <- function(RefBetasSub, CellTypes) {
  #transpose the refbeta sub matrix
  data = t(RefBetasSub)
  #the rows now correspond to the cell type (and subject) specific labels (samples)
  n = nrow(data)
  #the cols now represent the cpgs
  J = ncol(data)
  #get each uniqe cell type without the subject specific info
  batchlevs = levels(factor(CellTypes))
  #calculate the number of unique pairs of cell types there are
  numPairs = choose(length(batchlevs), 2)
  #get the pairs 
  Pairs = combn(batchlevs, 2)
  PairNames = apply(Pairs, 2, function(w) paste(w[1], " vs ", w[2], sep = ""))
  
  #calculate overall weights since the number of samples for each cell type
  #may differ depending on the data we have
  overallweights = table(factor(CellTypes))/length(CellTypes)
  weights = list()
  for(p in 1:(numPairs+1)) {
    if(p == 1) weights[[p]] = overallweights
    else {
      tmp = rep(0, length(batchlevs))
      tmp[batchlevs %in% Pairs[,(p-1)]] = 0.5
      names(tmp) = batchlevs
      weights[[p]] = tmp
    }
  }
  names(weights) = c("Overall", PairNames)
  
  #Now calulate the centroids for each cell type
  
  avg.w = NULL
  #number of cpgs by number of cell types
  cent.within = matrix(nrow = J, ncol = length(batchlevs))
  #name of the cps
  rownames(cent.within) = colnames(data)
  #name of the cell types
  colnames(cent.within) = batchlevs
  #for each cell type we will calculate the centriod
  for(i in 1:length(batchlevs)) {
    #get all of the columns that are of cell type i
    dat.w = as.matrix(data[CellTypes == batchlevs[i],])
    #if the number of columns of the matrix we just got has one column, then transpose
    #if there is only one sample then matrix algebra won't work unless transposed
    if (ncol(dat.w) == 1) dat.w = t(dat.w)
    #number of different rows we have of that cell type (samples)
    n.sub = nrow(dat.w)
    #get the averages across each of the n.sub rows (across the samples) for each cpg
    #should be a J x 1 vector
    cent.w = t(dat.w) %*% rep(1/n.sub, n.sub)
    #add this to our initialized matrix
    #add to all the rows of the ith column 
    cent.within[,i] = as.vector(cent.w)
    
    #for each sample of cell type i 
    #note dopar means do in parallel 
    dist.w = foreach(j=1:nrow(dat.w)) %do% {
      #let x be the jth row 
      x = as.numeric(dat.w[j,])
      #get a vector of differences between that sample and the average accross samples
      d = x-cent.w
      #now square those values, add them together and take the square root of that value
      sqrt(t(d) %*% d)
    }
    #get the average of these distances 
    avg.w[i] = mean(unlist(dist.w), na.rm = T)
  }
  #name of the cell types
  names(avg.w) = batchlevs
  
  #weights is a list, apply this to each element of that list
  DSCvals = lapply(weights, function(r) { 
    #each element of what weights contain (here weights for each cell type)
    wt = r
    #unique cell type where the wt not equal to 0
    cells = batchlevs[wt!=0]
    #get an index for each sample
    samp.ind = CellTypes %in% cells
    #get the data for each of the above indexes (for each sample)
    data.ind = data[samp.ind,]
    #get the number of rows
    n.ind = nrow(data.ind)
    #things will start to differ here
    #get the overall centroid 
    cent.global.ind = t(data.ind) %*% rep(1/n.ind, n.ind)
    #
    cent.within.ind = cent.within[,cells]
    #THIS IS WHAT WE WANT TO CHANGE
    dist.b = apply(cent.within.ind, 2, function(v) {
      dist.ind = v - as.vector(cent.global.ind)
      sqrt(t(dist.ind) %*% dist.ind)
    })
    #get the between cell type disperson
    Db = wt[cells] %*% dist.b
    #get the within cell type dispersion
    #average for within cluster dispersion for the cells that are being compared
    Dw = wt[cells] %*% avg.w[cells]
    DSC = Db/Dw
  })	
  
  DSCvals = unlist(DSCvals)	
  return(DSCvals)
}


####################################################################################
#FUNCTION:   distances.v1
####################################################################################
# ARGUMENTS:
#
# 	RefBetasSub:    A U x n matrix of methylation beta-values for samples 
#                   where DNAm was profiled in isolated leukocyte subtypes 
#                   U indicates the number of candidate DMRs and n
#                   indicates the number of samples 
#
# 	CellTypes:      A vector of length n whose elements are the cell 
#                   identities of the corresponding to the columns
#                   of RefBetasSub 
#########################################################################
# OUTPUT:        	Within cluster (cell type) distances
#                    
#########################################################################


distances.v1 <- function(RefBetasSub, CellTypes) {
  #transpose the refbeta sub matrix
  data = t(RefBetasSub)
  #the rows now correspond to the cell type (and subject) specific labels (samples)
  n = nrow(data)
  #the cols now represent the cpgs
  J = ncol(data)
  #get each uniqe cell type without the subject specific info
  batchlevs = levels(factor(CellTypes))
  #calculate the number of unique pairs of cell types there are
  numPairs = choose(length(batchlevs), 2)
  #get the pairs 
  Pairs = combn(batchlevs, 2)
  PairNames = apply(Pairs, 2, function(w) paste(w[1], " vs ", w[2], sep = ""))
  
  #calculate overall weights since the number of samples for each cell type
  #may differ depending on the data we have
  overallweights = table(factor(CellTypes))/length(CellTypes)
  
  #Now calulate the centroids for each cell type
  
  avg.w = NULL
  #number of cpgs by number of cell types
  cent.within = matrix(nrow = J, ncol = length(batchlevs))
  #name of the cps
  rownames(cent.within) = colnames(data)
  #name of the cell types
  colnames(cent.within) = batchlevs
  #for each cell type we will calculate the centriod
  for(i in 1:length(batchlevs)) {
    #get all of the columns that are of cell type i
    dat.w = as.matrix(data[CellTypes == batchlevs[i],])
    #if the number of columns of the matrix we just got has one column, then transpose
    #if there is only one sample then matrix algebra won't work unless transposed
    if (ncol(dat.w) == 1) dat.w = t(dat.w)
    #number of different rows we have of that cell type (samples)
    n.sub = nrow(dat.w)
    #get the averages across each of the n.sub rows (across the samples) for each cpg
    #should be a J x 1 vector
    cent.w = t(dat.w) %*% rep(1/n.sub, n.sub)
    #add this to our initialized matrix
    #add to all the rows of the ith column 
    cent.within[,i] = as.vector(cent.w)
    
    #for each sample of cell type i 
    #note dopar means do in parallel 
    dist.w = foreach(j=1:nrow(dat.w)) %do% {
      #let x be the jth row 
      x = as.numeric(dat.w[j,])
      #get a vector of differences between that sample and the average accross samples
      d = x-cent.w
      #now square those values, add them together and take the square root of that value
      sqrt(t(d) %*% d)
    }
    #get the average of these distances 
    avg.w[i] = mean(unlist(dist.w), na.rm = T)
  }
  #name of the cell types
  names(avg.w) = batchlevs
  
  
  
  #get the between cluster distances
  pair.dist = as.vector(dist(t(cent.within), method = 'euclidean'))
  names(pair.dist) = PairNames
  
  out=list(within = avg.w, pairwise = pair.dist)
  
  return(out)
  
}


#####################################################################################
# Function for updating probabilities of a CpG being selected (IN PARALLEL)
#####################################################################################
#FUNCTION ARGUMENTS
#
#cpg_list: candidate vector of cpgs to randomly select from
#
#probabilities: vector of probabilities of selecting the cpgs in cpg_list
#
#referencedata: A U x n matrix of methylation beta-values for samples where DNAm was profiled 
#in isolated leukocyte subtypes U indicates the number of candidate DMRs and n
#indicates the number of samples  
#
#fixed: A T/F indicating whether the number of cpgs to be selected should be fixed
#
#n: The number of cpgs to randomly select if fixed=T
#
#lower: The lower limit of cpgs to randomly select if fixed=F
#
#upper: The upper limit of cpgs to randomly select if fixed=F
#
#cellident: vector of column identitites that tells us the cell typ associated with each column of
#the reference beta matrix (will be used to calculate distances of interest)
#
#FUNCTION RETURNS
#
#1) List of randomly selected CpGs
#2) vector of updated probabilities of a CpG being selected
#3) The new metric generated by the library with the randomly selected CpGs
#####################################################################################
# parallelized function
#####################################################################################

update_probabilities_par <- function(cpg_list, p, referencedata, cellident, fixed, n, lower, upper, ncores=2){
  #registerDoParallel(ncores)
  
  #if the number of cpgs to be randomly selected is fixed
  if(fixed==T){
    fixedcpgs <- sample(cpg_list, size=n, replace = F, prob = p)
    #subset reference betas to the cpgs slected from above
    fixedref <- referencedata[fixedcpgs,]
    
    #get distances of interest
    fixed_dist <- distances.v1(fixedref, cellident)
    
    #calculate new metric with all cpgs
    fixed_metric <- min(unlist(fixed_dist[2]))/min(unlist(fixed_dist[1]))
    
    #now check value of new metric after removing each cpg and update probabilities
    newp <- foreach(i = 1:length(fixedcpgs), .packages = 'foreach', .export = 'distances.v1') %dopar% {
      #remove cpg i from our matrix of betas
      cpgswithout <- fixedcpgs[-i]
      without <- referencedata[cpgswithout,]
      
      #calculate distances of interest without cpg i
      dist_without <- distances.v1(without, cellident)
      
      #calculate metric without cpg i
      metric_without<- min(unlist(dist_without[2]))/min(unlist(dist_without[1]))
      
      #calculate the ratio of the metric without the ith cpg and with all cpgs
      #will be used to update probabilities
      r <- metric_without/fixed_metric
      
      #if metric is smaller without this cpg then increase prob of this cpg being selected
      #if metric is bigger without this cpg then decrease prob of this cpg being selected
      if(metric_without != fixed_metric){
        temp <- (1/r)*p[fixedcpgs[i]]
      }
      #else the probability stays the same
      if(metric_without == fixed_metric){
        temp <- p[fixedcpgs[i]]
      }
      temp
    }
    #return probabilities, library and metric
    out=list(library = fixedcpgs, probability = unlist(newp), metric=fixed_metric)
    return(out)
    
  }
  #if the number of cpgs to be randomly selected is random
  if(fixed==F){
    #first get the random number of cpgs to be selected using the upper and lower limits given
    randomnum <- sample(lower:upper, 1)
    randomcpgs <- sample(cpg_list, size=randomnum, replace = F, prob = p)
    #subset reference betas to the cpgs slected from above
    randomref <- referencedata[randomcpgs,]
    
    #get distances of interest
    random_dist <- distances.v1(randomref, cellident)
    
    #calculate new metric with all cpgs
    random_metric <- min(unlist(random_dist[2]))/min(unlist(random_dist[1]))
    
    #now check value of new metric after removing each cpg and update probabilities
    newp <- foreach(i = 1:length(randomcpgs), .packages = 'foreach', .export = 'distances.v1') %dopar% {
      #remove cpg i from our matrix of betas
      cpgswithout <- randomcpgs[-i]
      without <- referencedata[cpgswithout,]
      
      #calculate distances of interest without cpg i
      dist_without <- distances.v1(without, cellident)
      
      #calculate metric without cpg i
      metric_without<- min(unlist(dist_without[2]))/min(unlist(dist_without[1]))
      
      #calculate the ratio of the metric without the ith cpg and with all cpgs
      #will be used to update probabilities
      r <- metric_without/random_metric
      
      
      #if metric is smaller without this cpg then increase prob of this cpg being selected
      #if metric is bigger without this cpg then decrease prob of this cpg being selected
      if(metric_without != random_metric){
        temp <- (1/r)*p[randomcpgs[i]]
      }
      #else the probability stays the same
      if(metric_without == random_metric){
        temp <- p[randomcpgs[i]]
      }
      temp
    }
    #return probabilities, library and metric
    out=list(library = randomcpgs, probability = unlist(newp), metric=random_metric)
    return(out)
  }
  
}
