# Simulation study to assess the ramifications of errors in cell predictions on EWAS 
# Type 1 error rate
# Devin C. Koestler, Ph.D.
# April 29th, 2015

#setwd("/Users/dkoestler/Documents/Mixture methodology/Composite/")

# load necessary packages
library(betareg)
library(gtools)
library(mvtnorm)

source("multiRegMultistat.R")
load("ObjectsForSimulationFDR.RData")
load("PhiEstimates.RData")
load("SimulationErrorRates.RData")

# Define the logit2 (base 2) and inverse-logit (base 2) functions
logit2 = function(x) log2(x) - log2(1-x)
expit2 = function(x) 2^x/(1 + 2^x)
logit = function(x) log(x) - log(1-x)

# Define the simulation parameters
J = nrow(Mmatrix)				 # total number of arrayed CpGs
G = 1000                        # number of probes selected for estimating FDR
n1 = c(100, 500)    # sample size in group 1
n2 = c(100, 500)	 # sample size in group 2
EucDistRange = 101               # range of the euclidean distance computed between the
								 # underlying cell composition of group (1,2) samps

#MSEs from the modified DSC results
MSEs_DSC <- matrix(NA, nrow = 2, ncol = 6, byrow = T)
rownames(MSEs_DSC) <- c("Legacy","Modified DSC")
colnames(MSEs_DSC) <- c("CD8", "CD4", "NK", "Bcell", "Mono", "Gran")

DSCmse <- c(0.0015658987, 0.0003365367,  0.0014113489,
            0.0001391400, 0.0002655767, 0.0004448588)
legmse <- c(0.0026881254,0.0016200706, 0.0028424012, 
            0.0001746676, 0.0003898145, 0.0006199289)

MSEs_DSC[1,] <- legmse
MSEs_DSC[2,] <- DSCmse

########################################################################################
# Step 0a: randomly sample G probes and subset the reference/target methylation matrix
########################################################################################
probes = sample(rownames(Mmatrix), G, replace = T)
M.star = as.matrix(Mmatrix[probes,])
Y.star = as.matrix(Y.full[probes,])
pi.ests = phiEsts[probes]

ind = 1
ErrorRates = list()
for(i in 1:length(n1)) {
		FDRs = matrix(NA, nrow = EucDistRange, ncol = 4)
		rownames(FDRs) = rownames(Alpha2)
		colnames(FDRs) = c("Truth", "Legacy", "Modified DSC", "Nothing")
		for(k in 1:EucDistRange) {
			N1 = n1[i]
			N2 = n2[i]
			N = N1 + N2	
			group = c(rep(1, N1), rep(2, N2))
			alpha2 = Alpha2[k,]

########################################################################################
# Step 1: Randomly sample mixture proportions for groups 1 and 2
########################################################################################
			omega.1 = rdirichlet(N1, alpha1)
			omega.2 = rdirichlet(N2, alpha2)
			omega = rbind(omega.1, omega.2)						 


########################################################################################
# Step 2: create the linear predictors for each subject and across each G CpGs and 
#         simulate the DNAm data (return both beta and M-values)
########################################################################################

		simBetas = matrix(NA, nrow = N, ncol = G)
		simMvals = matrix(NA, nrow = N, ncol = G)
		for(p in 1:N) {
			alpha.mu = as.vector(M.star%*%omega[p,])
			alpha.phi = pi.ests
	
			# if(group[p] == 1) {
			# 	eta.mu = logit(alpha.mu)
			# 	eta.phi = alpha.phi
			# }
			# else {
			# 	eta.mu = logit(alpha.mu) 
			# 	eta.phi = alpha.phi
			# }

			mu.star = alpha.mu
			phi.star = alpha.phi

			a =  mu.star * phi.star
			b = (1-mu.star) * phi.star
			for(g in 1:G) {
				simBetas[p,g] = rbeta(1, shape1 = a[g], shape2 = b[g])
				simMvals[p,g] = logit2(simBetas[p,g])
			}
		}

########################################################################################
# Step 3: Randomly generate mixture predictions based on MSPE estimated using the
#         various L-DMR subsets, i.e., estimateCellCounts, ours, ...
########################################################################################
		numMods = 2
		Omegas = list()
		for(m in 1:numMods) {
			omega.pred = matrix(NA, nrow = N, ncol = 6)
			for(p in 1:N) {
				omega.pred[p,] = rmvnorm(1, mean = omega[p,], sigma = diag(MSEs_DSC[m,]))
			}	

			omega.pred[omega.pred <= 0] = 0
			omega.pred = t(apply(omega.pred, 1, function(w) {
					x.sum = sum(w);
					return(w/x.sum)}))
			Omegas[[m]] = omega.pred
		}
		names(Omegas) = c("Legacy", "Modified DSC")

########################################################################################
# Step 5: Do the EWAS to estimate type 1 error rate
########################################################################################

		Contrast = list()
		Contrast[[1]] = matrix(c(0,1,0,0,0,0,0), nrow = 7)
		Contrast1 = list()
		Contrast1[[1]] = matrix(c(0,1), nrow = 2)

		DesignMatTruth = model.matrix(~factor(group, levels = c(2,1)) + omega[,-6])
		DesignMatLegacy = model.matrix(~factor(group, levels = c(2,1)) + Omegas[[1]][,-6])
		DesignMatDSC = model.matrix(~factor(group, levels = c(2,1)) + Omegas[[2]][,-6])
		DesignMatNothing = model.matrix(~factor(group, levels = c(2,1)))

		fitTruth = multiRegMultistat(t(simMvals), DesignMatTruth, Llist = Contrast)
		fitLegacy = multiRegMultistat(t(simMvals), DesignMatLegacy, Llist = Contrast)
		fitDSC  = multiRegMultistat(t(simMvals), DesignMatDSC, Llist = Contrast)
		fitNothing = multiRegMultistat(t(simMvals), DesignMatNothing, Llist = Contrast1)

########################################################################################
# Step 6: Save the type 1 error rates
########################################################################################
		FDRs[k,] = c(sum(fitTruth$pvalue <= 0.05)/G,
		 			sum(fitLegacy$pvalue <= 0.05)/G,
		 			sum(fitDSC$pvalue <= 0.05)/G,
		 			sum(fitNothing$pvalue <= 0.05)/G)
		}
		print(paste("Iteration: ", ind, " Finished", sep = ""))
		ErrorRates[[ind]] = FDRs
		#names(ErrorRates[[ind]]) = paste("(n1 = ", N1, "; n2 = ", N2, ")", sep = "")
		ind = ind + 1
}


########################################################################################
#plot the results
########################################################################################

library(reshape2)

fdrs100 <- data.frame(ErrorRates[[1]])
fdrs100$dis <- seq(0,100)
fdrs500 <- data.frame(ErrorRates[[2]])
fdrs500$dis <- seq(0,100)

####################
# 100 fdrs
####################

fdrs100plot <- ggplot(data = fdrs100, aes(x=dis, group= 1)) + 
                      geom_smooth(aes(y=Nothing, colour= "No Adjustment"), method = "loess", se=FALSE, size=1.5)+
                      geom_smooth(aes(y=Legacy, colour= "Legacy"), method = "loess", se=FALSE, size=1.5) +
                      geom_smooth(aes(y=Modified.DSC, colour= "RESET"), method = "loess", se=FALSE, size=1.5) +
                      geom_smooth(aes(y=Truth, colour= "Truth"), method = "loess", se=FALSE, size=1.5)+
                      theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                            legend.title = element_blank(), legend.position = c(0.25, 0.8),
                            legend.key.size = unit(.6, 'cm'), legend.text = element_text(size=15),
                            axis.text=element_text(size=17))


####################
# 100 fdrs diffs
####################

diff100 <-  fdrs100$Legacy - fdrs100$Modified.DSC
diff100dat <- data.frame(diff100)
diff100dat$dis <- seq(0,100)

#mean diff
mean.diff100 <- mean(diff100)

#loess smoother
diff100plot <- ggplot(data = diff100dat, aes(x=dis, y=diff100)) + 
                geom_smooth(method = "loess", se=FALSE, size=1.5, color="black")+
              geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1.5)+
              geom_hline(yintercept=mean.diff100, linetype="dashed", color = "red", size=1.5)+
              theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                    axis.text=element_text(size=17))

#no loess smoother
# diff100plot2 <- ggplot(data = diff100dat, aes(x=dis, y=diff100)) + 
#   geom_line()+
#   geom_hline(yintercept=0, linetype="dashed", color = "black", size=1.5)+
#   geom_hline(yintercept=mean.diff100, linetype="dashed", color = "red", size=1.5)+
#   theme(axis.title.x = element_blank(), axis.title.y = element_blank())

####################
# 500 fdrs
####################

fdrs500plot <- ggplot(data = fdrs500, aes(x=dis, group= 1)) + 
  geom_smooth(aes(y=Nothing, colour= "No Adjustment"), method = "loess", se=FALSE, size=1.5)+
  geom_smooth(aes(y=Legacy, colour= "Legacy"), method = "loess", se=FALSE, size=1.5) +
  geom_smooth(aes(y=Modified.DSC, colour= "RESET"), method = "loess", se=FALSE, size=1.5) +
  geom_smooth(aes(y=Truth, colour= "Truth"), method = "loess", se=FALSE, size=1.5)+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.title = element_blank(), legend.position = c(0.25, 0.8),
        legend.key.size = unit(.6, 'cm'), legend.text = element_text(size=15),
        axis.text=element_text(size=17))

####################
# 500 fdrs diffs
####################

diff500 <-  fdrs500$Legacy - fdrs500$Modified.DSC
diff500dat <- data.frame(diff500)
diff500dat$dis <- seq(0,100)

#mean diff
mean.diff500 <- mean(diff500)

#loess smoother
diff500plot <- ggplot(data = diff500dat, aes(x=dis, y=diff500)) + 
  geom_smooth(method = "loess", se=FALSE, size=1.5, color="black")+
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1.5)+
  geom_hline(yintercept=mean.diff500, linetype="dashed", color = "red", size=1.5)+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text=element_text(size=17))

#no loess smoother
# diff500plot2 <- ggplot(data = diff500dat, aes(x=dis, y=diff500)) + 
#   geom_line()+
#   geom_hline(yintercept=0, linetype="dashed", color = "black", size=1.5)+
#   geom_hline(yintercept=mean.diff500, linetype="dashed", color = "red", size=1.5)+
#   theme(axis.title.x = element_blank(), axis.title.y = element_blank())


############################
# put together on one plot
############################

library(gridExtra)

grid.arrange(fdrs100plot, fdrs500plot, diff100plot,
             diff500plot, nrow=2)


