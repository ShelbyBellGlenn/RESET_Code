#load results from model fits and cell prop predictions
load("Cell Preds DSC Feinberg.RData")
load("Cell Preds DSC Hannum.RData")
load("Cell Preds Legacy Feinberg.RData")
load("Cell Preds Legacy Hannum.RData")
load("Feinberg Legacy Results.RData")
load("Feinberg Modified DSC Results.RData")
load("Hannum Legacy Results.RData")
load("Hannum modifed DSC Results.RData")

#packages
library(ggplot2)
library(ggrepel)
library(gridExtra)

##############################################################
#Feinberg figures 120
##############################################################

###############
#Histogram
###############

DSCFeinR2 = MultiDSCFein$R2
LegFeinR2 = MultiLegFein$R2
NothingR2Fein = rep(0, length(LegFeinR2))
R2DiffFein = DSCFeinR2 - LegFeinR2

hist(R2DiffFein, col = "grey", cex.axis = 1.8, ylim = c(0,140000))
mean(R2DiffFein>0)

fein <- data.frame(R2DiffFein)
feinhist <- ggplot(fein, aes(x=R2DiffFein))+ 
  geom_histogram(bins=20,color="Black", fill= "gray") +
  theme_classic() +
  labs(y="Frequency", x= expression("Difference in "~R^{2}))


######################################
#Scatter plot of estimates of leg DSC
######################################

cd4pred <- data.frame(CellPredsDSCFein[,"CD4T"], CellPredslegFein[,"CD4T"])
cd4plot <- ggplot(data = cd4pred, aes(x=CellPredslegFein[,"CD4T"], y=CellPredsDSCFein[,"CD4T"])) + 
  geom_point(color="blue", size=1) + 
  ggtitle("CD4") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                               axis.title.x = element_blank(), 
                               axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="blue")


cd8pred <- data.frame(CellPredsDSCFein[,"CD8T"], CellPredslegFein[,"CD8T"])
cd8plot <- ggplot(data = cd8pred, aes(x=CellPredslegFein[,"CD8T"], y=CellPredsDSCFein[,"CD8T"])) + 
  geom_point(color="red", size=1) + 
  ggtitle("CD8") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                               axis.title.x = element_blank(), 
                               axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="red")

bpred <- data.frame(CellPredsDSCFein[,"Bcell"], CellPredslegFein[,"Bcell"])
bplot <- ggplot(data = bpred, aes(x=CellPredslegFein[,"Bcell"], y=CellPredsDSCFein[,"Bcell"])) + 
  geom_point(color="dark green", size=1) + 
  ggtitle("B") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                               axis.title.x = element_blank(), 
                               axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="dark green")


nkpred <- data.frame(CellPredsDSCFein[,"NK"], CellPredslegFein[,"NK"])
nkplot <- ggplot(data = nkpred, aes(x=CellPredslegFein[,"NK"], y=CellPredsDSCFein[,"NK"])) + 
  geom_point(color="purple", size=1) + 
  ggtitle("NK") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                             axis.title.x = element_blank(), 
                             axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="purple")

monopred <- data.frame(CellPredsDSCFein[,"Mono"], CellPredslegFein[,"Mono"])
monoplot <- ggplot(data = monopred, aes(x=CellPredslegFein[,"Mono"], y=CellPredsDSCFein[,"Mono"])) + 
  geom_point(color="orange", size=1) + 
  ggtitle("Monocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                             axis.title.x = element_blank(), 
                             axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="orange")

granpred <- data.frame(CellPredsDSCFein[,"Gran"], CellPredslegFein[,"Gran"])
granplot <- ggplot(data = granpred, aes(x=CellPredslegFein[,"Gran"], y=CellPredsDSCFein[,"Gran"])) + 
  geom_point(color="deeppink", size=1) + 
  ggtitle("Granulocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                    axis.title.x = element_blank(), 
                                    axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="deeppink")


#put all 6 plots on one plot
grid.arrange(cd4plot, cd8plot, bplot,
             nkplot, monoplot, granplot, nrow=2)



##############################################################
#Hannum figures 120
##############################################################

DSCHanR2 = MultiDSCHan$R2
LegHanR2 = MultilegHan$R2
NothingR2han = rep(0, length(LegHanR2))
R2DiffHan = DSCHanR2 - LegHanR2


hist(R2DiffHan, col = "grey", cex.axis = 1.8, ylim = c(0,140000))
mean(R2DiffHan>0)

han <- data.frame(R2DiffHan)
hanhist <- ggplot(han, aes(x=R2DiffHan))+ 
  geom_histogram(bins=20,color="Black", fill= "gray") +
  theme_classic() +
  labs(y="Frequency", x= expression("Difference in "~R^{2}))

######################################
#Scatter plot of estimates of leg DSC
######################################

cd4pred2 <- data.frame(CellPredsDSCHan[,"CD4T"], CellPredslegHan[,"CD4T"])
cd4plot2 <- ggplot(data = cd4pred2, aes(x=CellPredslegHan[,"CD4T"], y=CellPredsDSCHan[,"CD4T"])) + 
  geom_point(color="blue", size=1) + 
  ggtitle("CD4") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                               axis.title.x = element_blank(), 
                               axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="blue")


cd8pred2 <- data.frame(CellPredsDSCHan[,"CD8T"], CellPredslegHan[,"CD8T"])
cd8plot2 <- ggplot(data = cd8pred2, aes(x=CellPredslegHan[,"CD8T"], y=CellPredsDSCHan[,"CD8T"])) + 
  geom_point(color="red", size=1) + 
  ggtitle("CD8") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                               axis.title.x = element_blank(), 
                               axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="red")

bpred2 <- data.frame(CellPredsDSCHan[,"Bcell"], CellPredslegHan[,"Bcell"])
bplot2 <- ggplot(data = bpred2, aes(x=CellPredslegHan[,"Bcell"], y=CellPredsDSCHan[,"Bcell"])) + 
  geom_point(color="dark green", size=1) + 
  ggtitle("B") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                             axis.title.x = element_blank(), 
                             axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="dark green")


nkpred2 <- data.frame(CellPredsDSCHan[,"NK"], CellPredslegHan[,"NK"])
nkplot2 <- ggplot(data = nkpred2, aes(x=CellPredslegHan[,"NK"], y=CellPredsDSCHan[,"NK"])) + 
  geom_point(color="purple", size=1) + 
  ggtitle("NK") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="purple")

monopred2 <- data.frame(CellPredsDSCHan[,"Mono"], CellPredslegHan[,"Mono"])
monoplot2 <- ggplot(data = monopred2, aes(x=CellPredslegHan[,"Mono"], y=CellPredsDSCHan[,"Mono"])) + 
  geom_point(color="orange", size=1) + 
  ggtitle("Monocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                    axis.title.x = element_blank(), 
                                    axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="orange")

granpred2 <- data.frame(CellPredsDSCHan[,"Gran"], CellPredslegHan[,"Gran"])
granplot2 <- ggplot(data = granpred2, aes(x=CellPredslegHan[,"Gran"], y=CellPredsDSCHan[,"Gran"])) + 
  geom_point(color="deeppink", size=1) + 
  ggtitle("Granulocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                       axis.title.x = element_blank(), 
                                       axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="deeppink")


#put all 6 plots on one plot
grid.arrange(cd4plot2, cd8plot2, bplot2,
             nkplot2, monoplot2, granplot2, nrow=2)



##############################################################
#Feinberg figures 72
##############################################################

#load data 
load("Cell Preds DSC Feinberg 72.RData")
load("Cell Preds Legacy Feinberg 72.RData")
load("Feinberg Legacy Results 72.RData")
load("Feinberg Modified DSC Results 72.RData")

DSCFeinR272 = MultiDSCFein72$R2
LegFeinR272 = MultiLegFein72$R2
NothingR2Fein72 = rep(0, length(LegFeinR272))
R2DiffFein72 = DSCFeinR272 - LegFeinR272

mean(R2DiffFein72>0)

fein72 <- data.frame(R2DiffFein72)
feinhist72 <- ggplot(fein72, aes(x=R2DiffFein72))+ 
  geom_histogram(bins=20,color="Black", fill= "gray") +
  theme_classic() +
  labs(y="Frequency", x= expression("Difference in "~R^{2}))


cd4pred72 <- data.frame(CellPredsDSCFein72[,"CD4T"], CellPredslegFein72[,"CD4T"])
cd4plot72 <- ggplot(data = cd4pred72, aes(x=CellPredslegFein72[,"CD4T"], y=CellPredsDSCFein72[,"CD4T"])) + 
  geom_point(color="blue", size=1) + 
  ggtitle("CD4") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="blue")


cd8pred72 <- data.frame(CellPredsDSCFein72[,"CD8T"], CellPredslegFein72[,"CD8T"])
cd8plot72 <- ggplot(data = cd8pred72, aes(x=CellPredslegFein72[,"CD8T"], y=CellPredsDSCFein72[,"CD8T"])) + 
  geom_point(color="red", size=1) + 
  ggtitle("CD8") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="red")

bpred72 <- data.frame(CellPredsDSCFein72[,"Bcell"], CellPredslegFein72[,"Bcell"])
bplot72 <- ggplot(data = bpred72, aes(x=CellPredslegFein72[,"Bcell"], y=CellPredsDSCFein72[,"Bcell"])) + 
  geom_point(color="dark green", size=1) + 
  ggtitle("B") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                       axis.title.x = element_blank(), 
                       axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="dark green")


nkpred72 <- data.frame(CellPredsDSCFein72[,"NK"], CellPredslegFein72[,"NK"])
nkplot72 <- ggplot(data = nkpred72, aes(x=CellPredslegFein72[,"NK"], y=CellPredsDSCFein72[,"NK"])) + 
  geom_point(color="purple", size=1) + 
  ggtitle("NK") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                        axis.title.x = element_blank(), 
                        axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="purple")

monopred72 <- data.frame(CellPredsDSCFein72[,"Mono"], CellPredslegFein72[,"Mono"])
monoplot72 <- ggplot(data = monopred72, aes(x=CellPredslegFein72[,"Mono"], y=CellPredsDSCFein72[,"Mono"])) + 
  geom_point(color="orange", size=1) + 
  ggtitle("Monocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="orange")

granpred72 <- data.frame(CellPredsDSCFein72[,"Gran"], CellPredslegFein72[,"Gran"])
granplot72 <- ggplot(data = granpred72, aes(x=CellPredslegFein72[,"Gran"], y=CellPredsDSCFein72[,"Gran"])) + 
  geom_point(color="deeppink", size=1) + 
  ggtitle("Granulocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                 axis.title.x = element_blank(), 
                                 axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="deeppink")


#put all 6 plots on one plot
grid.arrange(cd4plot72, cd8plot72, bplot72,
             nkplot72, monoplot72, granplot72, nrow=2)


##############################################################
#Hannum figures 72
##############################################################

#load data
load("Cell Preds DSC Hannum 72.RData")
load("Cell Preds Legacy Hannum 72.RData")
load("Hannum Legacy Results 72.RData")
load("Hannum modifed DSC Results 72.RData")

DSCHanR272 = MultiDSCHan72$R2
LegHanR272 = MultilegHan72$R2
NothingR2han72 = rep(0, length(LegHanR272))
R2DiffHan72 = DSCHanR272 - LegHanR272

mean(R2DiffHan72>0)

han72 <- data.frame(R2DiffHan72)
hanhist72 <- ggplot(han72, aes(x=R2DiffHan72))+ 
  geom_histogram(bins=20,color="Black", fill= "gray") +
  theme_classic() +
  labs(y="Frequency", x= expression("Difference in "~R^{2}))

cd4pred2_72 <- data.frame(CellPredsDSCHan72[,"CD4T"], CellPredslegHan72[,"CD4T"])
cd4plot2_72 <- ggplot(data = cd4pred2_72, aes(x=CellPredslegHan72[,"CD4T"], y=CellPredsDSCHan72[,"CD4T"])) + 
  geom_point(color="blue", size=1) + 
  ggtitle("CD4") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="blue")


cd8pred2_72 <- data.frame(CellPredsDSCHan72[,"CD8T"], CellPredslegHan72[,"CD8T"])
cd8plot2_72 <- ggplot(data = cd8pred2_72, aes(x=CellPredslegHan72[,"CD8T"], y=CellPredsDSCHan72[,"CD8T"])) + 
  geom_point(color="red", size=1) + 
  ggtitle("CD8") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="red")

bpred2_72 <- data.frame(CellPredsDSCHan72[,"Bcell"], CellPredslegHan72[,"Bcell"])
bplot2_72 <- ggplot(data = bpred2_72, aes(x=CellPredslegHan72[,"Bcell"], y=CellPredsDSCHan72[,"Bcell"])) + 
  geom_point(color="dark green", size=1) + 
  ggtitle("B") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                       axis.title.x = element_blank(), 
                       axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="dark green")


nkpred2_72 <- data.frame(CellPredsDSCHan72[,"NK"], CellPredslegHan72[,"NK"])
nkplot2_72 <- ggplot(data = nkpred2_72, aes(x=CellPredslegHan72[,"NK"], y=CellPredsDSCHan72[,"NK"])) + 
  geom_point(color="purple", size=1) + 
  ggtitle("NK") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                        axis.title.x = element_blank(), 
                        axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="purple")

monopred2_72 <- data.frame(CellPredsDSCHan72[,"Mono"], CellPredslegHan72[,"Mono"])
monoplot2_72 <- ggplot(data = monopred2_72, aes(x=CellPredslegHan72[,"Mono"], y=CellPredsDSCHan72[,"Mono"])) + 
  geom_point(color="orange", size=1) + 
  ggtitle("Monocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="orange")

granpred2_72 <- data.frame(CellPredsDSCHan72[,"Gran"], CellPredslegHan72[,"Gran"])
granplot2_72 <- ggplot(data = granpred2_72, aes(x=CellPredslegHan72[,"Gran"], y=CellPredsDSCHan72[,"Gran"])) + 
  geom_point(color="deeppink", size=1) + 
  ggtitle("Granulocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                 axis.title.x = element_blank(), 
                                 axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="deeppink")


#put all 6 plots on one plot
grid.arrange(cd4plot2_72, cd8plot2_72, bplot2_72,
             nkplot2_72, monoplot2_72, granplot2_72, nrow=2)

##############################################################
#Feinberg figures 180
##############################################################

#load data 
load("Cell Preds DSC Feinberg 180.RData")
load("Cell Preds Legacy Feinberg 180.RData")
load("Feinberg Legacy Results 180.RData")
load("Feinberg Modified DSC Results 180.RData")

DSCFeinR2180 = MultiDSCFein180$R2
LegFeinR2180 = MultiLegFein180$R2
NothingR2Fein180 = rep(0, length(LegFeinR2180))
R2DiffFein180 = DSCFeinR2180 - LegFeinR2180

mean(R2DiffFein180>0)

fein180 <- data.frame(R2DiffFein180)
feinhist180 <- ggplot(fein180, aes(x=R2DiffFein180))+ 
  geom_histogram(bins=20,color="Black", fill= "gray") +
  theme_classic() +
  labs(y="Frequency", x= expression("Difference in "~R^{2}))

cd4pred180 <- data.frame(CellPredsDSCFein180[,"CD4T"], CellPredslegFein180[,"CD4T"])
cd4plot180 <- ggplot(data = cd4pred180, aes(x=CellPredslegFein180[,"CD4T"], y=CellPredsDSCFein180[,"CD4T"])) + 
  geom_point(color="blue", size=1) + 
  ggtitle("CD4") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="blue")


cd8pred180 <- data.frame(CellPredsDSCFein180[,"CD8T"], CellPredslegFein180[,"CD8T"])
cd8plot180 <- ggplot(data = cd8pred180, aes(x=CellPredslegFein180[,"CD8T"], y=CellPredsDSCFein180[,"CD8T"])) + 
  geom_point(color="red", size=1) + 
  ggtitle("CD8") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="red")

bpred180 <- data.frame(CellPredsDSCFein180[,"Bcell"], CellPredslegFein180[,"Bcell"])
bplot180 <- ggplot(data = bpred180, aes(x=CellPredslegFein180[,"Bcell"], y=CellPredsDSCFein180[,"Bcell"])) + 
  geom_point(color="dark green", size=1) + 
  ggtitle("B") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                       axis.title.x = element_blank(), 
                       axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="dark green")


nkpred180 <- data.frame(CellPredsDSCFein180[,"NK"], CellPredslegFein180[,"NK"])
nkplot180 <- ggplot(data = nkpred180, aes(x=CellPredslegFein180[,"NK"], y=CellPredsDSCFein180[,"NK"])) + 
  geom_point(color="purple", size=1) + 
  ggtitle("NK") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                        axis.title.x = element_blank(), 
                        axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="purple")

monopred180 <- data.frame(CellPredsDSCFein180[,"Mono"], CellPredslegFein180[,"Mono"])
monoplot180 <- ggplot(data = monopred180, aes(x=CellPredslegFein180[,"Mono"], y=CellPredsDSCFein180[,"Mono"])) + 
  geom_point(color="orange", size=1) + 
  ggtitle("Monocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="orange")

granpred180 <- data.frame(CellPredsDSCFein180[,"Gran"], CellPredslegFein180[,"Gran"])
granplot180 <- ggplot(data = granpred180, aes(x=CellPredslegFein180[,"Gran"], y=CellPredsDSCFein180[,"Gran"])) + 
  geom_point(color="deeppink", size=1) + 
  ggtitle("Granulocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                 axis.title.x = element_blank(), 
                                 axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="deeppink")


#put all 6 plots on one plot
grid.arrange(cd4plot180, cd8plot180, bplot180,
             nkplot180, monoplot180, granplot180, nrow=2)

##############################################################
#Hannum figures 180
##############################################################

#load data
load("Cell Preds DSC Hannum 180.RData")
load("Cell Preds Legacy Hannum 180.RData")
load("Hannum Legacy Results 180.RData")
load("Hannum modifed DSC Results 180.RData")

DSCHanR2180 = MultiDSCHan180$R2
LegHanR2180 = MultilegHan180$R2
NothingR2han180 = rep(0, length(LegHanR2180))
R2DiffHan180 = DSCHanR2180 - LegHanR2180

mean(R2DiffHan180>0)

han180 <- data.frame(R2DiffHan180)
hanhist180 <- ggplot(han180, aes(x=R2DiffHan180))+ 
  geom_histogram(bins=20,color="Black", fill= "gray") +
  theme_classic() +
  labs(y="Frequency", x= expression("Difference in "~R^{2}))

cd4pred2_180 <- data.frame(CellPredsDSCHan180[,"CD4T"], CellPredslegHan180[,"CD4T"])
cd4plot2_180 <- ggplot(data = cd4pred2_180, aes(x=CellPredslegHan180[,"CD4T"], y=CellPredsDSCHan180[,"CD4T"])) + 
  geom_point(color="blue", size=1) + 
  ggtitle("CD4") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="blue")


cd8pred2_180 <- data.frame(CellPredsDSCHan180[,"CD8T"], CellPredslegHan180[,"CD8T"])
cd8plot2_180 <- ggplot(data = cd8pred2_180, aes(x=CellPredslegHan180[,"CD8T"], y=CellPredsDSCHan180[,"CD8T"])) + 
  geom_point(color="red", size=1) + 
  ggtitle("CD8") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="red")

bpred2_180 <- data.frame(CellPredsDSCHan180[,"Bcell"], CellPredslegHan180[,"Bcell"])
bplot2_180 <- ggplot(data = bpred2_180, aes(x=CellPredslegHan180[,"Bcell"], y=CellPredsDSCHan180[,"Bcell"])) + 
  geom_point(color="dark green", size=1) + 
  ggtitle("B") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                       axis.title.x = element_blank(), 
                       axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="dark green")


nkpred2_180 <- data.frame(CellPredsDSCHan180[,"NK"], CellPredslegHan180[,"NK"])
nkplot2_180 <- ggplot(data = nkpred2_180, aes(x=CellPredslegHan180[,"NK"], y=CellPredsDSCHan180[,"NK"])) + 
  geom_point(color="purple", size=1) + 
  ggtitle("NK") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                        axis.title.x = element_blank(), 
                        axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="purple")

monopred2_180 <- data.frame(CellPredsDSCHan180[,"Mono"], CellPredslegHan180[,"Mono"])
monoplot2_180 <- ggplot(data = monopred2_180, aes(x=CellPredslegHan180[,"Mono"], y=CellPredsDSCHan180[,"Mono"])) + 
  geom_point(color="orange", size=1) + 
  ggtitle("Monocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="orange")

granpred2_180 <- data.frame(CellPredsDSCHan180[,"Gran"], CellPredslegHan180[,"Gran"])
granplot2_180 <- ggplot(data = granpred2_180, aes(x=CellPredslegHan180[,"Gran"], y=CellPredsDSCHan180[,"Gran"])) + 
  geom_point(color="deeppink", size=1) + 
  ggtitle("Granulocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                 axis.title.x = element_blank(), 
                                 axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="deeppink")


#put all 6 plots on one plot
grid.arrange(cd4plot2_180, cd8plot2_180, bplot2_180,
             nkplot2_180, monoplot2_180, granplot2_180, nrow=2)

##############################################################
#Feinberg figures 240
##############################################################

#load data 
load("Cell Preds DSC Feinberg 240.RData")
load("Cell Preds Legacy Feinberg 240.RData")
load("Feinberg Legacy Results 240.RData")
load("Feinberg Modified DSC Results 240.RData")

DSCFeinR2240 = MultiDSCFein240$R2
LegFeinR2240 = MultiLegFein240$R2
NothingR2Fein240 = rep(0, length(LegFeinR2240))
R2DiffFein240 = DSCFeinR2240 - LegFeinR2240

mean(R2DiffFein240>0)

fein240 <- data.frame(R2DiffFein240)
feinhist240 <- ggplot(fein240, aes(x=R2DiffFein240))+ 
  geom_histogram(bins=20,color="Black", fill= "gray") +
  theme_classic() +
  labs(y="Frequency", x= expression("Difference in "~R^{2}))

cd4pred240 <- data.frame(CellPredsDSCFein240[,"CD4T"], CellPredslegFein240[,"CD4T"])
cd4plot240 <- ggplot(data = cd4pred240, aes(x=CellPredslegFein240[,"CD4T"], y=CellPredsDSCFein240[,"CD4T"])) + 
  geom_point(color="blue", size=1) + 
  ggtitle("CD4") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="blue")


cd8pred240 <- data.frame(CellPredsDSCFein240[,"CD8T"], CellPredslegFein240[,"CD8T"])
cd8plot240 <- ggplot(data = cd8pred240, aes(x=CellPredslegFein240[,"CD8T"], y=CellPredsDSCFein240[,"CD8T"])) + 
  geom_point(color="red", size=1) + 
  ggtitle("CD8") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="red")

bpred240 <- data.frame(CellPredsDSCFein240[,"Bcell"], CellPredslegFein240[,"Bcell"])
bplot240 <- ggplot(data = bpred240, aes(x=CellPredslegFein240[,"Bcell"], y=CellPredsDSCFein240[,"Bcell"])) + 
  geom_point(color="dark green", size=1) + 
  ggtitle("B") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                       axis.title.x = element_blank(), 
                       axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="dark green")


nkpred240 <- data.frame(CellPredsDSCFein240[,"NK"], CellPredslegFein240[,"NK"])
nkplot240 <- ggplot(data = nkpred240, aes(x=CellPredslegFein240[,"NK"], y=CellPredsDSCFein240[,"NK"])) + 
  geom_point(color="purple", size=1) + 
  ggtitle("NK") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                        axis.title.x = element_blank(), 
                        axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="purple")

monopred240 <- data.frame(CellPredsDSCFein240[,"Mono"], CellPredslegFein240[,"Mono"])
monoplot240 <- ggplot(data = monopred240, aes(x=CellPredslegFein240[,"Mono"], y=CellPredsDSCFein240[,"Mono"])) + 
  geom_point(color="orange", size=1) + 
  ggtitle("Monocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="orange")

granpred240 <- data.frame(CellPredsDSCFein240[,"Gran"], CellPredslegFein240[,"Gran"])
granplot240 <- ggplot(data = granpred240, aes(x=CellPredslegFein240[,"Gran"], y=CellPredsDSCFein240[,"Gran"])) + 
  geom_point(color="deeppink", size=1) + 
  ggtitle("Granulocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                 axis.title.x = element_blank(), 
                                 axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="deeppink")


#put all 6 plots on one plot
grid.arrange(cd4plot240, cd8plot240, bplot240,
             nkplot240, monoplot240, granplot240, nrow=2)

##############################################################
#Hannum figures 240
##############################################################

#load data
load("Cell Preds DSC Hannum 240.RData")
load("Cell Preds Legacy Hannum 240.RData")
load("Hannum Legacy Results 240.RData")
load("Hannum modifed DSC Results 240.RData")

DSCHanR2240 = MultiDSCHan240$R2
LegHanR2240 = MultilegHan240$R2
NothingR2han240 = rep(0, length(LegHanR2240))
R2DiffHan240 = DSCHanR2240 - LegHanR2240

mean(R2DiffHan240>0)

han240 <- data.frame(R2DiffHan240)
hanhist240 <- ggplot(han240, aes(x=R2DiffHan240))+ 
  geom_histogram(bins=20,color="Black", fill= "gray") +
  theme_classic() +
  labs(y="Frequency", x= expression("Difference in "~R^{2}))

cd4pred2_240 <- data.frame(CellPredsDSCHan240[,"CD4T"], CellPredslegHan240[,"CD4T"])
cd4plot2_240 <- ggplot(data = cd4pred2_240, aes(x=CellPredslegHan240[,"CD4T"], y=CellPredsDSCHan240[,"CD4T"])) + 
  geom_point(color="blue", size=1) + 
  ggtitle("CD4") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="blue")


cd8pred2_240 <- data.frame(CellPredsDSCHan240[,"CD8T"], CellPredslegHan240[,"CD8T"])
cd8plot2_240 <- ggplot(data = cd8pred2_240, aes(x=CellPredslegHan240[,"CD8T"], y=CellPredsDSCHan240[,"CD8T"])) + 
  geom_point(color="red", size=1) + 
  ggtitle("CD8") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="red")

bpred2_240 <- data.frame(CellPredsDSCHan240[,"Bcell"], CellPredslegHan240[,"Bcell"])
bplot2_240 <- ggplot(data = bpred2_240, aes(x=CellPredslegHan240[,"Bcell"], y=CellPredsDSCHan240[,"Bcell"])) + 
  geom_point(color="dark green", size=1) + 
  ggtitle("B") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                       axis.title.x = element_blank(), 
                       axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="dark green")


nkpred2_240 <- data.frame(CellPredsDSCHan240[,"NK"], CellPredslegHan240[,"NK"])
nkplot2_240 <- ggplot(data = nkpred2_240, aes(x=CellPredslegHan240[,"NK"], y=CellPredsDSCHan240[,"NK"])) + 
  geom_point(color="purple", size=1) + 
  ggtitle("NK") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                        axis.title.x = element_blank(), 
                        axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="purple")

monopred2_240 <- data.frame(CellPredsDSCHan240[,"Mono"], CellPredslegHan240[,"Mono"])
monoplot2_240 <- ggplot(data = monopred2_240, aes(x=CellPredslegHan240[,"Mono"], y=CellPredsDSCHan240[,"Mono"])) + 
  geom_point(color="orange", size=1) + 
  ggtitle("Monocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="orange")

granpred2_240 <- data.frame(CellPredsDSCHan240[,"Gran"], CellPredslegHan240[,"Gran"])
granplot2_240 <- ggplot(data = granpred2_240, aes(x=CellPredslegHan240[,"Gran"], y=CellPredsDSCHan240[,"Gran"])) + 
  geom_point(color="deeppink", size=1) + 
  ggtitle("Granulocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                 axis.title.x = element_blank(), 
                                 axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="deeppink")


#put all 6 plots on one plot
grid.arrange(cd4plot2_240, cd8plot2_240, bplot2_240,
             nkplot2_240, monoplot2_240, granplot2_240, nrow=2)

##############################################################
#Feinberg figures 300
##############################################################

#load data 
load("Cell Preds DSC Feinberg 300.RData")
load("Cell Preds Legacy Feinberg 300.RData")
load("Feinberg Legacy Results 300.RData")
load("Feinberg Modified DSC Results 300.RData")

DSCFeinR2300 = MultiDSCFein300$R2
LegFeinR2300 = MultiLegFein300$R2
NothingR2Fein300 = rep(0, length(LegFeinR2300))
R2DiffFein300 = DSCFeinR2300 - LegFeinR2300

mean(R2DiffFein300>0)

fein300 <- data.frame(R2DiffFein300)
feinhist300 <- ggplot(fein300, aes(x=R2DiffFein300))+ 
  geom_histogram(bins=20,color="Black", fill= "gray") +
  theme_classic() +
  labs(y="Frequency", x= expression("Difference in "~R^{2}))

cd4pred300 <- data.frame(CellPredsDSCFein300[,"CD4T"], CellPredslegFein300[,"CD4T"])
cd4plot300 <- ggplot(data = cd4pred300, aes(x=CellPredslegFein300[,"CD4T"], y=CellPredsDSCFein300[,"CD4T"])) + 
  geom_point(color="blue", size=1) + 
  ggtitle("CD4") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="blue")


cd8pred300 <- data.frame(CellPredsDSCFein300[,"CD8T"], CellPredslegFein300[,"CD8T"])
cd8plot300 <- ggplot(data = cd8pred300, aes(x=CellPredslegFein300[,"CD8T"], y=CellPredsDSCFein300[,"CD8T"])) + 
  geom_point(color="red", size=1) + 
  ggtitle("CD8") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="red")

bpred300 <- data.frame(CellPredsDSCFein300[,"Bcell"], CellPredslegFein300[,"Bcell"])
bplot300 <- ggplot(data = bpred300, aes(x=CellPredslegFein300[,"Bcell"], y=CellPredsDSCFein300[,"Bcell"])) + 
  geom_point(color="dark green", size=1) + 
  ggtitle("B") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                       axis.title.x = element_blank(), 
                       axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="dark green")


nkpred300 <- data.frame(CellPredsDSCFein300[,"NK"], CellPredslegFein300[,"NK"])
nkplot300 <- ggplot(data = nkpred300, aes(x=CellPredslegFein300[,"NK"], y=CellPredsDSCFein300[,"NK"])) + 
  geom_point(color="purple", size=1) + 
  ggtitle("NK") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                        axis.title.x = element_blank(), 
                        axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="purple")

monopred300 <- data.frame(CellPredsDSCFein300[,"Mono"], CellPredslegFein300[,"Mono"])
monoplot300 <- ggplot(data = monopred300, aes(x=CellPredslegFein300[,"Mono"], y=CellPredsDSCFein300[,"Mono"])) + 
  geom_point(color="orange", size=1) + 
  ggtitle("Monocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="orange")

granpred300 <- data.frame(CellPredsDSCFein300[,"Gran"], CellPredslegFein300[,"Gran"])
granplot300 <- ggplot(data = granpred300, aes(x=CellPredslegFein300[,"Gran"], y=CellPredsDSCFein300[,"Gran"])) + 
  geom_point(color="deeppink", size=1) + 
  ggtitle("Granulocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                 axis.title.x = element_blank(), 
                                 axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="deeppink")


#put all 6 plots on one plot
grid.arrange(cd4plot300, cd8plot300, bplot300,
             nkplot300, monoplot300, granplot300, nrow=2)

##############################################################
#Hannum figures 300
##############################################################

#load data
load("Cell Preds DSC Hannum 300.RData")
load("Cell Preds Legacy Hannum 300.RData")
load("Hannum Legacy Results 300.RData")
load("Hannum modifed DSC Results 300.RData")

DSCHanR2300 = MultiDSCHan300$R2
LegHanR2300 = MultilegHan300$R2
NothingR2han300 = rep(0, length(LegHanR2300))
R2DiffHan300 = DSCHanR2300 - LegHanR2300

mean(R2DiffHan300>0)

han300 <- data.frame(R2DiffHan300)
hanhist300 <- ggplot(han300, aes(x=R2DiffHan300))+ 
  geom_histogram(bins=20,color="Black", fill= "gray") +
  theme_classic() +
  labs(y="Frequency", x= expression("Difference in "~R^{2}))

cd4pred2_300 <- data.frame(CellPredsDSCHan300[,"CD4T"], CellPredslegHan300[,"CD4T"])
cd4plot2_300 <- ggplot(data = cd4pred2_300, aes(x=CellPredslegHan300[,"CD4T"], y=CellPredsDSCHan300[,"CD4T"])) + 
  geom_point(color="blue", size=1) + 
  ggtitle("CD4") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="blue")


cd8pred2_300 <- data.frame(CellPredsDSCHan300[,"CD8T"], CellPredslegHan300[,"CD8T"])
cd8plot2_300 <- ggplot(data = cd8pred2_300, aes(x=CellPredslegHan300[,"CD8T"], y=CellPredsDSCHan300[,"CD8T"])) + 
  geom_point(color="red", size=1) + 
  ggtitle("CD8") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="red")

bpred2_300 <- data.frame(CellPredsDSCHan300[,"Bcell"], CellPredslegHan300[,"Bcell"])
bplot2_300 <- ggplot(data = bpred2_300, aes(x=CellPredslegHan300[,"Bcell"], y=CellPredsDSCHan300[,"Bcell"])) + 
  geom_point(color="dark green", size=1) + 
  ggtitle("B") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                       axis.title.x = element_blank(), 
                       axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="dark green")


nkpred2_300 <- data.frame(CellPredsDSCHan300[,"NK"], CellPredslegHan300[,"NK"])
nkplot2_300 <- ggplot(data = nkpred2_300, aes(x=CellPredslegHan300[,"NK"], y=CellPredsDSCHan300[,"NK"])) + 
  geom_point(color="purple", size=1) + 
  ggtitle("NK") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                        axis.title.x = element_blank(), 
                        axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="purple")

monopred2_300 <- data.frame(CellPredsDSCHan300[,"Mono"], CellPredslegHan300[,"Mono"])
monoplot2_300 <- ggplot(data = monopred2_300, aes(x=CellPredslegHan300[,"Mono"], y=CellPredsDSCHan300[,"Mono"])) + 
  geom_point(color="orange", size=1) + 
  ggtitle("Monocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="orange")

granpred2_300 <- data.frame(CellPredsDSCHan300[,"Gran"], CellPredslegHan300[,"Gran"])
granplot2_300 <- ggplot(data = granpred2_300, aes(x=CellPredslegHan300[,"Gran"], y=CellPredsDSCHan300[,"Gran"])) + 
  geom_point(color="deeppink", size=1) + 
  ggtitle("Granulocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                 axis.title.x = element_blank(), 
                                 axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="deeppink")


#put all 6 plots on one plot
grid.arrange(cd4plot2_300, cd8plot2_300, bplot2_300,
             nkplot2_300, monoplot2_300, granplot2_300, nrow=2)

##############################################################
#Feinberg figures 360
##############################################################

#load data 
load("Cell Preds DSC Feinberg 360.RData")
load("Cell Preds Legacy Feinberg 360.RData")
load("Feinberg Legacy Results 360.RData")
load("Feinberg Modified DSC Results 360.RData")

DSCFeinR2360 = MultiDSCFein360$R2
LegFeinR2360 = MultiLegFein360$R2
NothingR2Fein360 = rep(0, length(LegFeinR2360))
R2DiffFein360 = DSCFeinR2360 - LegFeinR2360

mean(R2DiffFein360>0)

fein360 <- data.frame(R2DiffFein360)
feinhist360 <- ggplot(fein360, aes(x=R2DiffFein360))+ 
  geom_histogram(bins=20,color="Black", fill= "gray") +
  theme_classic() +
  labs(y="Frequency", x= expression("Difference in "~R^{2}))

cd4pred360 <- data.frame(CellPredsDSCFein360[,"CD4T"], CellPredslegFein360[,"CD4T"])
cd4plot360 <- ggplot(data = cd4pred360, aes(x=CellPredslegFein360[,"CD4T"], y=CellPredsDSCFein360[,"CD4T"])) + 
  geom_point(color="blue", size=1) + 
  ggtitle("CD4") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="blue")


cd8pred360 <- data.frame(CellPredsDSCFein360[,"CD8T"], CellPredslegFein360[,"CD8T"])
cd8plot360 <- ggplot(data = cd8pred360, aes(x=CellPredslegFein360[,"CD8T"], y=CellPredsDSCFein360[,"CD8T"])) + 
  geom_point(color="red", size=1) + 
  ggtitle("CD8") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="red")

bpred360 <- data.frame(CellPredsDSCFein360[,"Bcell"], CellPredslegFein360[,"Bcell"])
bplot360 <- ggplot(data = bpred360, aes(x=CellPredslegFein360[,"Bcell"], y=CellPredsDSCFein360[,"Bcell"])) + 
  geom_point(color="dark green", size=1) + 
  ggtitle("B") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                       axis.title.x = element_blank(), 
                       axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="dark green")


nkpred360 <- data.frame(CellPredsDSCFein360[,"NK"], CellPredslegFein360[,"NK"])
nkplot360 <- ggplot(data = nkpred360, aes(x=CellPredslegFein360[,"NK"], y=CellPredsDSCFein360[,"NK"])) + 
  geom_point(color="purple", size=1) + 
  ggtitle("NK") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                        axis.title.x = element_blank(), 
                        axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="purple")

monopred360 <- data.frame(CellPredsDSCFein360[,"Mono"], CellPredslegFein360[,"Mono"])
monoplot360 <- ggplot(data = monopred360, aes(x=CellPredslegFein360[,"Mono"], y=CellPredsDSCFein360[,"Mono"])) + 
  geom_point(color="orange", size=1) + 
  ggtitle("Monocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="orange")

granpred360 <- data.frame(CellPredsDSCFein360[,"Gran"], CellPredslegFein360[,"Gran"])
granplot360 <- ggplot(data = granpred360, aes(x=CellPredslegFein360[,"Gran"], y=CellPredsDSCFein360[,"Gran"])) + 
  geom_point(color="deeppink", size=1) + 
  ggtitle("Granulocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                 axis.title.x = element_blank(), 
                                 axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="deeppink")


#put all 6 plots on one plot
grid.arrange(cd4plot360, cd8plot360, bplot360,
             nkplot360, monoplot360, granplot360, nrow=2)

##############################################################
#Hannum figures 360
##############################################################

#load data
load("Cell Preds DSC Hannum 360.RData")
load("Cell Preds Legacy Hannum 360.RData")
load("Hannum Legacy Results 360.RData")
load("Hannum modifed DSC Results 360.RData")

DSCHanR2360 = MultiDSCHan360$R2
LegHanR2360 = MultilegHan360$R2
NothingR2han360 = rep(0, length(LegHanR2360))
R2DiffHan360 = DSCHanR2360 - LegHanR2360

mean(R2DiffHan360>0)

han360 <- data.frame(R2DiffHan360)
hanhist360 <- ggplot(han360, aes(x=R2DiffHan360))+ 
  geom_histogram(bins=20,color="Black", fill= "gray") +
  theme_classic() +
  labs(y="Frequency", x= expression("Difference in "~R^{2}))

cd4pred2_360 <- data.frame(CellPredsDSCHan360[,"CD4T"], CellPredslegHan360[,"CD4T"])
cd4plot2_360 <- ggplot(data = cd4pred2_360, aes(x=CellPredslegHan360[,"CD4T"], y=CellPredsDSCHan360[,"CD4T"])) + 
  geom_point(color="blue", size=1) + 
  ggtitle("CD4") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="blue")


cd8pred2_360 <- data.frame(CellPredsDSCHan360[,"CD8T"], CellPredslegHan360[,"CD8T"])
cd8plot2_360 <- ggplot(data = cd8pred2_360, aes(x=CellPredslegHan360[,"CD8T"], y=CellPredsDSCHan360[,"CD8T"])) + 
  geom_point(color="red", size=1) + 
  ggtitle("CD8") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="red")

bpred2_360 <- data.frame(CellPredsDSCHan360[,"Bcell"], CellPredslegHan360[,"Bcell"])
bplot2_360 <- ggplot(data = bpred2_360, aes(x=CellPredslegHan360[,"Bcell"], y=CellPredsDSCHan360[,"Bcell"])) + 
  geom_point(color="dark green", size=1) + 
  ggtitle("B") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                       axis.title.x = element_blank(), 
                       axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="dark green")


nkpred2_360 <- data.frame(CellPredsDSCHan360[,"NK"], CellPredslegHan360[,"NK"])
nkplot2_360 <- ggplot(data = nkpred2_360, aes(x=CellPredslegHan360[,"NK"], y=CellPredsDSCHan360[,"NK"])) + 
  geom_point(color="purple", size=1) + 
  ggtitle("NK") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                        axis.title.x = element_blank(), 
                        axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="purple")

monopred2_360 <- data.frame(CellPredsDSCHan360[,"Mono"], CellPredslegHan360[,"Mono"])
monoplot2_360 <- ggplot(data = monopred2_360, aes(x=CellPredslegHan360[,"Mono"], y=CellPredsDSCHan360[,"Mono"])) + 
  geom_point(color="orange", size=1) + 
  ggtitle("Monocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="orange")

granpred2_360 <- data.frame(CellPredsDSCHan360[,"Gran"], CellPredslegHan360[,"Gran"])
granplot2_360 <- ggplot(data = granpred2_360, aes(x=CellPredslegHan360[,"Gran"], y=CellPredsDSCHan360[,"Gran"])) + 
  geom_point(color="deeppink", size=1) + 
  ggtitle("Granulocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                 axis.title.x = element_blank(), 
                                 axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="deeppink")


#put all 6 plots on one plot
grid.arrange(cd4plot2_360, cd8plot2_360, bplot2_360,
             nkplot2_360, monoplot2_360, granplot2_360, nrow=2)

##############################################################
#Feinberg figures 540
##############################################################

#load data 
load("Cell Preds DSC Feinberg 540.RData")
load("Cell Preds Legacy Feinberg 540.RData")
load("Feinberg Legacy Results 540.RData")
load("Feinberg Modified DSC Results 540.RData")

DSCFeinR2540 = MultiDSCFein540$R2
LegFeinR2540 = MultiLegFein540$R2
NothingR2Fein540 = rep(0, length(LegFeinR2540))
R2DiffFein540 = DSCFeinR2540 - LegFeinR2540

mean(R2DiffFein540>0)

fein540 <- data.frame(R2DiffFein540)
feinhist540 <- ggplot(fein540, aes(x=R2DiffFein540))+ 
  geom_histogram(bins=20,color="Black", fill= "gray") +
  theme_classic() +
  labs(y="Frequency", x= expression("Difference in "~R^{2}))

cd4pred540 <- data.frame(CellPredsDSCFein540[,"CD4T"], CellPredslegFein540[,"CD4T"])
cd4plot540 <- ggplot(data = cd4pred540, aes(x=CellPredslegFein540[,"CD4T"], y=CellPredsDSCFein540[,"CD4T"])) + 
  geom_point(color="blue", size=1) + 
  ggtitle("CD4") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="blue")


cd8pred540 <- data.frame(CellPredsDSCFein540[,"CD8T"], CellPredslegFein540[,"CD8T"])
cd8plot540 <- ggplot(data = cd8pred540, aes(x=CellPredslegFein540[,"CD8T"], y=CellPredsDSCFein540[,"CD8T"])) + 
  geom_point(color="red", size=1) + 
  ggtitle("CD8") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="red")

bpred540 <- data.frame(CellPredsDSCFein540[,"Bcell"], CellPredslegFein540[,"Bcell"])
bplot540 <- ggplot(data = bpred540, aes(x=CellPredslegFein540[,"Bcell"], y=CellPredsDSCFein540[,"Bcell"])) + 
  geom_point(color="dark green", size=1) + 
  ggtitle("B") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                       axis.title.x = element_blank(), 
                       axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="dark green")


nkpred540 <- data.frame(CellPredsDSCFein540[,"NK"], CellPredslegFein540[,"NK"])
nkplot540 <- ggplot(data = nkpred540, aes(x=CellPredslegFein540[,"NK"], y=CellPredsDSCFein540[,"NK"])) + 
  geom_point(color="purple", size=1) + 
  ggtitle("NK") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                        axis.title.x = element_blank(), 
                        axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="purple")

monopred540 <- data.frame(CellPredsDSCFein540[,"Mono"], CellPredslegFein540[,"Mono"])
monoplot540 <- ggplot(data = monopred540, aes(x=CellPredslegFein540[,"Mono"], y=CellPredsDSCFein540[,"Mono"])) + 
  geom_point(color="orange", size=1) + 
  ggtitle("Monocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="orange")

granpred540 <- data.frame(CellPredsDSCFein540[,"Gran"], CellPredslegFein540[,"Gran"])
granplot540 <- ggplot(data = granpred540, aes(x=CellPredslegFein540[,"Gran"], y=CellPredsDSCFein540[,"Gran"])) + 
  geom_point(color="deeppink", size=1) + 
  ggtitle("Granulocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                 axis.title.x = element_blank(), 
                                 axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="deeppink")


#put all 6 plots on one plot
grid.arrange(cd4plot540, cd8plot540, bplot540,
             nkplot540, monoplot540, granplot540, nrow=2)

##############################################################
#Hannum figures 540
##############################################################

#load data
load("Cell Preds DSC Hannum 540.RData")
load("Cell Preds Legacy Hannum 540.RData")
load("Hannum Legacy Results 540.RData")
load("Hannum modifed DSC Results 540.RData")

DSCHanR2540 = MultiDSCHan540$R2
LegHanR2540 = MultilegHan540$R2
NothingR2han540 = rep(0, length(LegHanR2540))
R2DiffHan540 = DSCHanR2540 - LegHanR2540

mean(R2DiffHan540>0)

han540 <- data.frame(R2DiffHan540)
hanhist540 <- ggplot(han540, aes(x=R2DiffHan540))+ 
  geom_histogram(bins=20,color="Black", fill= "gray") +
  theme_classic() +
  labs(y="Frequency", x= expression("Difference in "~R^{2}))

cd4pred2_540 <- data.frame(CellPredsDSCHan540[,"CD4T"], CellPredslegHan540[,"CD4T"])
cd4plot2_540 <- ggplot(data = cd4pred2_540, aes(x=CellPredslegHan540[,"CD4T"], y=CellPredsDSCHan540[,"CD4T"])) + 
  geom_point(color="blue", size=1) + 
  ggtitle("CD4") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="blue")


cd8pred2_540 <- data.frame(CellPredsDSCHan540[,"CD8T"], CellPredslegHan540[,"CD8T"])
cd8plot2_540 <- ggplot(data = cd8pred2_540, aes(x=CellPredslegHan540[,"CD8T"], y=CellPredsDSCHan540[,"CD8T"])) + 
  geom_point(color="red", size=1) + 
  ggtitle("CD8") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="red")

bpred2_540 <- data.frame(CellPredsDSCHan540[,"Bcell"], CellPredslegHan540[,"Bcell"])
bplot2_540 <- ggplot(data = bpred2_540, aes(x=CellPredslegHan540[,"Bcell"], y=CellPredsDSCHan540[,"Bcell"])) + 
  geom_point(color="dark green", size=1) + 
  ggtitle("B") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                       axis.title.x = element_blank(), 
                       axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="dark green")


nkpred2_540 <- data.frame(CellPredsDSCHan540[,"NK"], CellPredslegHan540[,"NK"])
nkplot2_540 <- ggplot(data = nkpred2_540, aes(x=CellPredslegHan540[,"NK"], y=CellPredsDSCHan540[,"NK"])) + 
  geom_point(color="purple", size=1) + 
  ggtitle("NK") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                        axis.title.x = element_blank(), 
                        axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="purple")

monopred2_540 <- data.frame(CellPredsDSCHan540[,"Mono"], CellPredslegHan540[,"Mono"])
monoplot2_540 <- ggplot(data = monopred2_540, aes(x=CellPredslegHan540[,"Mono"], y=CellPredsDSCHan540[,"Mono"])) + 
  geom_point(color="orange", size=1) + 
  ggtitle("Monocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="orange")

granpred2_540 <- data.frame(CellPredsDSCHan540[,"Gran"], CellPredslegHan540[,"Gran"])
granplot2_540 <- ggplot(data = granpred2_540, aes(x=CellPredslegHan540[,"Gran"], y=CellPredsDSCHan540[,"Gran"])) + 
  geom_point(color="deeppink", size=1) + 
  ggtitle("Granulocyte") + theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                 axis.title.x = element_blank(), 
                                 axis.title.y = element_blank()) +
  geom_abline(size=1) + geom_smooth(method=lm, se=FALSE, color="deeppink")


#put all 6 plots on one plot
grid.arrange(cd4plot2_540, cd8plot2_540, bplot2_540,
             nkplot2_540, monoplot2_540, granplot2_540, nrow=2)

