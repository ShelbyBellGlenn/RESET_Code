###########################################################################
# Proof of principle for paper
##########################################################################

#load needed functions
source("UnsupervisedOptimalDMRFinderV3.R")

# column names for data (3 cell types, 6 samples)
samples <- c("Cell One_1","Cell One_2","Cell One_3",
             "Cell One_4","Cell One_5","Cell One_6",
             "Cell Two_1","Cell Two_2","Cell Two_3",
             "Cell Two_4","Cell Two_5","Cell Two_6",
             "Cell Three_1","Cell Three_2","Cell Three_3",
             "Cell Three_4","Cell Three_5","Cell Three_6")

#get the cell type identities
cell_ident <- sapply(strsplit(samples, split="_", fixed = T), function(x) (x[1]))

#number of features we would like to simulate data for (e.g. cpgs)
num_features <- 100

#percentage of features we would like to exhibit a difference
perc_features <- seq(from=0, to=100, by=5)


#name of the features
features <- paste0("Feature ", seq(1,100))

#magnitude of the difference we would like to exhibit
magnitude <- c(1:10)

#create empty matrix to store the results: rows will be the % of features which
#exhibit any difference and the columns will be the magnitude of that difference
DSC_mat <- matrix(NA, nrow = length(perc_features), ncol = length(magnitude))
rownames(DSC_mat) <- paste0(seq(from=0, to=100, by=5), "%")
colnames(DSC_mat) <- seq(1,10)

#store the data genereated at each iteration of filling the DSC matrix
data.matrices <- NULL
#index to keep track of where the next data matrix should go
index <- 1

#note the "base" normal distribution will be a normal(0,1)
for(i in 1:length(perc_features)){
  for(j in 1:length(magnitude)){
    #set up data matrix to fill at each iteration of this loop
    data.mat <- matrix(NA, nrow = length(features), ncol = length(samples))
    rownames(data.mat) <- features
    colnames(data.mat) <- samples
    
    if(i>1){
      #generate data for the number of features which will exhbit a difference
      for(k in 1:perc_features[i]){
        #randomly select which cell type will have which mean for features which exhibit a difference
        delta <- c(0, magnitude[j], -1*magnitude[j])
        means <- sample(delta, 3, replace = F)
        
        #means[1] = mean for the features which will exhibit a difference for cell type 1
        
        #generate data for cell type 1
        data.mat[k,1] <- rnorm(1, means[1], 1)
        data.mat[k,2] <- rnorm(1, means[1], 1)
        data.mat[k,3] <- rnorm(1, means[1], 1)
        data.mat[k,4] <- rnorm(1, means[1], 1)
        data.mat[k,5] <- rnorm(1, means[1], 1)
        data.mat[k,6] <- rnorm(1, means[1], 1)
        
        #generate data for cell type 2
        data.mat[k,7] <- rnorm(1, means[2], 1)
        data.mat[k,8] <- rnorm(1, means[2], 1)
        data.mat[k,9] <- rnorm(1, means[2], 1)
        data.mat[k,10] <- rnorm(1, means[2], 1)
        data.mat[k,11] <- rnorm(1, means[2], 1)
        data.mat[k,12] <- rnorm(1, means[2], 1)
        
        #generate data for cell type 3
        data.mat[k,13] <- rnorm(1, means[3], 1)
        data.mat[k,14] <- rnorm(1, means[3], 1)
        data.mat[k,15] <- rnorm(1, means[3], 1)
        data.mat[k,16] <- rnorm(1, means[3], 1)
        data.mat[k,17] <- rnorm(1, means[3], 1)
        data.mat[k,18] <- rnorm(1, means[3], 1)
      }
    }
    
    #genereate data for the features which will not exhibit a difference
    #generate data for cell type 1
    if(i < 21){
      num <- perc_features[i] + 1
      data.mat[num:100,1] <- rnorm(100-perc_features[i], 0, 1)
      data.mat[num:100,2] <- rnorm(100-perc_features[i], 0, 1) 
      data.mat[num:100,3] <- rnorm(100-perc_features[i], 0, 1) 
      data.mat[num:100,4] <- rnorm(100-perc_features[i], 0, 1) 
      data.mat[num:100,5] <- rnorm(100-perc_features[i], 0, 1) 
      data.mat[num:100,6] <- rnorm(100-perc_features[i], 0, 1) 
      
      #generate data for cell type 2
      data.mat[num:100,7] <- rnorm(100-perc_features[i], 0, 1) 
      data.mat[num:100,8] <- rnorm(100-perc_features[i], 0, 1) 
      data.mat[num:100,9] <- rnorm(100-perc_features[i], 0, 1) 
      data.mat[num:100,10] <- rnorm(100-perc_features[i], 0, 1) 
      data.mat[num:100,11] <- rnorm(100-perc_features[i], 0, 1) 
      data.mat[num:100,12] <- rnorm(100-perc_features[i], 0, 1) 
      
      #generate data for cell type 3
      data.mat[num:100,13] <- rnorm(100-perc_features[i], 0, 1) 
      data.mat[num:100,14] <- rnorm(100-perc_features[i], 0, 1) 
      data.mat[num:100,15] <- rnorm(100-perc_features[i], 0, 1) 
      data.mat[num:100,16] <- rnorm(100-perc_features[i], 0, 1) 
      data.mat[num:100,17] <- rnorm(100-perc_features[i], 0, 1) 
      data.mat[num:100,18] <- rnorm(100-perc_features[i], 0, 1) 
    }
    
    #store the randomly generated data matrix for this iteration
    data.matrices[[index]] <- data.mat
    index <- index + 1
    
    #calcualte the DSC
    dists <- distances.v1(data.mat, cell_ident)
    DSC_check <- min(unlist(dists[2]))/min(unlist(dists[1]))
    
    #store that value in our matrix of DSC values
    DSC_mat[i,j] <- DSC_check
  }
}

#visualize results
library(gplots)
library(RColorBrewer)

m <- c(5,5)
heatmap.2(DSC_mat, col = colorRampPalette(c("pink", "red"))(32), 
          trace = "n", dendrogram = "n", Rowv = F, Colv = F, margins = m, 
          xlab = "Magnitude of Difference Between Features",
          ylab="Percent of Features Exhibiting Differences",
          density.info = "none", key.xlab = "DSC Value", key.title = "none" ,
          lmat=rbind(c(4, 2), c(1, 3)), lhei=c(2, 8), lwid=c(4, 1))  


#####################################################
# create heatmaps for three scenarios 
#####################################################

#want colors for columns corresponding to cell type
cell.cols = brewer.pal(3,"Pastel1")[c(3,2,1)]       
sample.cols = c(rep(cell.cols[1], 6),
                rep(cell.cols[2], 6),
                rep(cell.cols[3], 6) )  

m2 <- c(7,5)
heatmap.2(data.matrices[[1]], Rowv = T, Colv = F, dendrogram = "row", 
          col = colorRampPalette(c("yellow", "black", "blue"))(64), trace = "none", 
          ColSideColors = sample.cols, margins=m2, main = "0% of Features Different with a 
          Magnitude of 1")

heatmap.2(data.matrices[[105]], Rowv = T, Colv = F, dendrogram = "row", 
          col = colorRampPalette(c("yellow", "black", "blue"))(64), trace = "none", 
          ColSideColors = sample.cols, margins=m2, main = "50% of Features Different with a 
          Magnitude of 5")


heatmap.2(data.matrices[[210]], Rowv = T, Colv = F, dendrogram = "row", 
          col = colorRampPalette(c("yellow", "black", "blue"))(64), trace = "none", 
          ColSideColors = sample.cols, margins=m2, main= "100% of Features Different with a 
          Magnitude of 10")


####################################################################
#PCA
####################################################################

#load packages
library(ggplot2)
library(grid)
library(gridExtra)

#no differences in none of the features
pca.dat1 <- as.data.frame(t(data.matrices[[1]]))
pca1 <- prcomp(pca.dat1)

pca1_out <- as.data.frame(pca1$x)
pca1_out$group <- c("Cell Type One","Cell Type One","Cell Type One","Cell Type One","Cell Type One",
                    "Cell Type One", "Cell Type Two", "Cell Type Two","Cell Type Two",
                    "Cell Type Two","Cell Type Two","Cell Type Two","Cell Type Three",
                    "Cell Type Three","Cell Type Three","Cell Type Three","Cell Type Three","Cell Type Three")


pca.plot1 <-ggplot(pca1_out,aes(x=PC1,y=PC2,color=group )) +geom_point(size=4.25) + 
  theme(legend.position = c(.85, .83), legend.key.size = unit(1.5, 'cm'), 
        legend.title=element_blank(), legend.text = element_text(size = 15.5))+
  scale_color_manual(values = c("Cell Type One" = "#7CAE00", 
                                "Cell Type Two" = "#00BFC4", "Cell Type Three" = "#F8766D"))+
  theme(axis.text=element_text(size=17), axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17))

#50% differences in half of the features
pca.dat2 <- as.data.frame(t(data.matrices[[105]]))
pca2 <- prcomp(pca.dat2)

pca2_out <- as.data.frame(pca2$x)
pca2_out$group <- c("Cell Type One","Cell Type One","Cell Type One","Cell Type One","Cell Type One",
                    "Cell Type One", "Cell Type Two", "Cell Type Two","Cell Type Two",
                    "Cell Type Two","Cell Type Two","Cell Type Two","Cell Type Three",
                    "Cell Type Three","Cell Type Three","Cell Type Three","Cell Type Three","Cell Type Three")

pca.plot2 <-ggplot(pca2_out,aes(x=PC1,y=PC2,color=group ))+geom_point(size=4.25) + 
  theme(legend.position = c(.85, .83),legend.key.size = unit(1.5, 'cm'),
        legend.title=element_blank(), legend.text = element_text(size = 15.5)) +
  scale_color_manual(values = c("Cell Type One" = "#7CAE00", 
                                "Cell Type Two" = "#00BFC4", "Cell Type Three" = "#F8766D"))+
  theme(axis.text=element_text(size=17), axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17))

#100% differences in all of the features
pca.dat3 <- as.data.frame(t(data.matrices[[210]]))
pca3 <- prcomp(pca.dat3)

pca3_out <- as.data.frame(pca3$x)
pca3_out$group <- c("Cell Type One","Cell Type One","Cell Type One","Cell Type One","Cell Type One",
                    "Cell Type One", "Cell Type Two", "Cell Type Two","Cell Type Two",
                    "Cell Type Two","Cell Type Two","Cell Type Two","Cell Type Three",
                    "Cell Type Three","Cell Type Three","Cell Type Three","Cell Type Three","Cell Type Three")

pca.plot3<-ggplot(pca3_out,aes(x=PC1,y=PC2,color=group )) +
  geom_point(size=4.25) + theme(legend.position = c(.85, .83), legend.key.size = unit(1.5, 'cm'),
            legend.title=element_blank(), legend.text = element_text(size = 15.5)) +
  scale_color_manual(values = c("Cell Type One" = "#7CAE00", 
                                "Cell Type Two" = "#00BFC4", "Cell Type Three" = "#F8766D"))+
  theme(axis.text=element_text(size=17), axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17))

#########################################################################################################
#bar plots for cell composition
#########################################################################################################

library(dplyr)

load("MethodA.RData")
load("MethodB.RData")
load("AdultMixed.RData")

#method a

types <- rep( c("CD4", "CD8", "Bcell", "NK", "Monocyte", "Granulocyte"), 6)
samples <- c(rep("1", 6), rep("2", 6), rep("3", 6), rep("4", 6), rep("5", 6), rep("6", 6) )
prop <- c(13, 11, 16, 12, 23, 25,
          7, 19, 19, 15, 19, 21,
          6, 33, 8, 11, 19, 23,
          16,29, 7, 15, 22, 11,
          11,20, 20, 22, 10, 17,
          18,13, 26, 15, 22, 6)

methoda <- data.frame(samples, types, prop)

#method b

prop <- c(13, 2, 1, 4, 5, 75,
           16, 11, 1, 2, 7, 63,
           9, 6, 2, 0, 10, 73,
           14, 8, 2, 3, 6, 67,
           12, 5, 6, 7, 4, 66,
           15, 4, 4, 2, 5, 70)
methodb <- data.frame(samples, types, prop)

#adult mixed

prop <- c(24.37, 12.31, 6.62, 2.46, 5.49, 39.96,
           11.96, 5.87, 4.26, 2.75, 6.04, 66.19,
           18.1, 9.86, 4.78, 4.75, 5.21, 46.57,
           16.19, 15.26, 5.91, 4.81, 8.67, 44.06,
           18.17, 6.17, 1.77, 1.94, 6.14, 59.06,
           11, 4.72, 3.96, 2.62, 5.96, 67.83)
adultmixed <- data.frame(samples, types, prop)


#put into one plot together
data_together <- rbind(methoda, methodb, adultmixed)
data_together$data_set <- c(rep("Method A",36), rep("Method B",36), rep("Adult Mixed",36))

prop_toether <- ggplot(data = data_together, aes(fill=types, y=prop, x=samples))+ 
  geom_bar(position="stack", stat="identity" ) + 
  scale_fill_manual("Cell Types", values = c("CD4" = "blue", 
                                             "CD8" = "red2", "Bcell" = "dark green",
                                             "NK" = "purple", "Monocyte" = "orange",
                                             "Granulocyte"= "deeppink")) +
  labs(y = "Proportion", x= "Samples") + facet_wrap(~data_set) +
  theme(axis.text=element_text(size=17), strip.text.x = element_text(size = 16),
        legend.key.size = unit(1.5, 'cm'), legend.text = element_text(size = 14),
        legend.title = element_text(size=14), axis.title.x = element_blank(),
        axis.title.y = element_blank())
  


################################################################################################
# box plots by cell type and the two methods
################################################################################################

#put data together
library.size <- c(rep("72", 6), rep("120", 6), rep("180", 6), 
                  rep("240", 6), rep("300", 6),
                  rep("360", 6), rep("540", 6),
                  rep("72", 6), rep("120", 6), rep("180", 6), 
                  rep("240", 6), rep("300", 6),
                  rep("360", 6), rep("540", 6))
Method <- c(rep("RESET", 42), rep("Legacy", 42))
rsq <-c(0.867, 0.941, 0.784, 0.970, 0.952, 0.992,
        0.889, 0.951, 0.781, 0.980, 0.965, 0.992,
        0.907, 0.957, 0.823, 0.989, 0.948, 0.993,
        0.928, 0.948, 0.828, 0.992, 0.959, 0.995,
        0.911, 0.950, 0.838, 0.991, 0.968, 0.996,
        0.918, 0.948, 0.838, 0.994, 0.966, 0.995,
        0.921,  0.934, 0.891, 0.993, 0.960, 0.994,
        
        0.843, 0.829, 0.923, 0.976, 0.866, 0.993, 
        0.850, 0.838, 0.919, 0.984, 0.888, 0.996,
        0.846, 0.882, 0.930, 0.991, 0.914, 0.997,
        0.892, 0.885, 0.934, 0.993, 0.938, 0.996,
        0.888, 0.892, 0.925, 0.994, 0.945, 0.996,
        0.892, 0.894, 0.915, 0.994, 0.945, 0.966,
        0.933, 0.926, 0.880, 0.994, 0.954, 0.966)
rmse <- c(2.500, 3.936, 4.419, 1.361, 1.552, 2.860,
          2.221, 3.668, 3.980, 2.033, 1.290, 2.997,     
          2.148, 3.567, 3.539, 1.294, 1.991, 2.329,
          1.748, 3.513, 3.527, 0.843, 1.921, 1.981,
          1.807, 3.200, 3.639, 0.746, 1.553, 1.941,
          1.828, 3.627, 3.375, 0.895, 1.481, 2.320,
          2.108, 4.806, 4.33, 0.702, 1.705, 1.797,
          
          3.537, 4.904, 6.471, 1.951, 2.498, 2.813,
          3.965, 6.317, 6.847, 1.923, 2.432, 1.950,
          5.902, 5.402, 4.839, 0.941, 2.218, 1.879,
          4.668, 5.301, 4.854, 0.887, 1.712, 1.669,
          4.418, 4.759, 4.482, 0.910, 1.605, 1.705,
          4.068, 4.705, 4.944, 0.890, 1.657, 1.733,
          3.003, 4.066, 4.349, 0.943, 1.439, 1.545)


box.data <- data.frame(rsq, rmse, library.size, Method)
box.data$library.size <- factor(box.data$library.size,
                         c("72", "120", "180", "240",
                           "300", "360", "540"))

rsq.plot <- ggplot(box.data, aes(x=library.size, y=rsq, fill=Method)) + 
  geom_boxplot() + labs(y = expression(R^{2}), x= "Library Size") +
  stat_summary(fun=mean, geom="point", aes(group=Method), position=position_dodge(.8),
               color="black", size=2) + theme(axis.text=element_text(size=17),
                                              axis.title.x = element_blank(),
                                              axis.title.y = element_blank(),
                                              legend.key.size = unit(1.5, 'cm'),
                                              legend.text = element_text(size = 15),
                                              legend.title = element_text(size=15))

rmse.plot <- ggplot(box.data, aes(x=library.size, y=rmse, fill=Method)) + 
  geom_boxplot() + labs(y = "RMSE", x= "Library Size") +
  stat_summary(fun=mean, geom="point", aes(group=Method), position=position_dodge(.8),
               color="black", size=2) + theme(axis.text=element_text(size=18),
                                              axis.title.x = element_blank(),
                                              axis.title.y = element_blank(),
                                              legend.key.size = unit(1.5, 'cm'),
                                              legend.text = element_text(size = 17),
                                              legend.title = element_text(size=17))


