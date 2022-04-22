library(ggplot2)
library(factoextra)
library(FactoMineR)
library(gridExtra)
library(tidyverse)
library(cluster)
library(RColorBrewer)
library(plotly)

#Adding New Ages to my dataframe

Chubdata <- read_csv("ChubLake_fulldata_Jan26_clrXRF.csv")
age.depth <- read_csv("Chub_22_new_164_ageDepth.csv")
ChubMarch <- merge(Chubdata,age.depth,by="composite_depth_cm.x")

#XRF Data
#This space plots XRF elemental count data against age and creates grids of multiple plots 


XRF <- ChubMarch %>% select(Age_median,
                            Ca, Fe, Mn, Ti, K, Si, S, Al)

XRF <- XRF %>% drop_na()

Ca <- XRF %>% ggplot(aes(Ca, composite_depth_cm.x))+geom_path()+geom_point(pch =21, size = 2,  fill = "seagreen") + 
  scale_y_reverse(breaks = seq(0,10500,500)) + theme_bw() + ylab("Age kBP") + geom_hline(yintercept = c(9712,8200,4250,2800))

Si<- XRF %>% ggplot(aes(Si, Age_median))+geom_path()+geom_point(pch =21, size = 2,  fill = "seagreen")+ scale_y_reverse(breaks = seq(0,10500,500)) + theme_bw() + ylab("")+ geom_hline(yintercept = c(9712,8200,4250,2800))

Ti <- XRF %>% ggplot(aes(Ti, Age_median))+geom_path()+geom_point(pch =21, size = 2,  fill = "seagreen") + 
  scale_y_reverse(breaks = seq(0,10500,500)) + theme_bw() +  ylab("")+ geom_hline(yintercept = c(9712,8200,4250,2800))

K <- XRF %>% ggplot(aes(K, Age_median))+geom_path()+geom_point(pch =21, size = 2,  fill = "seagreen") + 
  scale_y_reverse(breaks = seq(0,10500,500)) + theme_bw() +  ylab("")+ geom_hline(yintercept = c(9712,8200,4250,2800))

grid.arrange(Ca, Si, Ti, K, nrow=1)

#PCA 
# This space clusters data using K-means clustering, and runs principal component analysis
# PCA is only run on XRF data 

PCA_ <- ChubMarch %>% select(Age_median, Al, Si, K, Ti, Fe, Mn, S, Ca, Sr, Zr, IncCoh)
PCA_ <- na.omit(PCA_)

# k-means clustering
# We use 3 k-clusters for this project, this section of code establishes the clusters
PCA_distance <- get_dist(PCA_[,2:11])
fviz_dist(PCA_distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

k3 <- kmeans(PCA_[,2:11], centers = 3, nstart = 25)

p2 <- fviz_cluster(k3, geom = "point",  data = PCA_[,2:11]) + ggtitle("k = 3")


set.seed(123)

## A quick test of plotting Si v Ca using kclusters as a factor, ensuring they ran correctly 
PCA_ <- PCA_ %>% mutate(
  kclusters = k3$cluster
)

PCA_ %>% ggplot(aes(S, Ca, fill = as.factor(kclusters))) + 
  geom_point(pch =21, size = 3) + 
  theme_bw() + 
  scale_fill_brewer(palette = "Greens") 

## re run to add k-clusters to PCA result 
pca.result <- prcomp(PCA_[,2:11])

fviz_eig(pca.result)
fviz_pca_ind(pca.result)
fviz_pca_ind(pca.result, geom.ind = "point", pointshape = 21, 
             pointsize = 3, 
             fill.ind = "grey", 
             col.ind = "black", 
             palette = "jco", 
             addEllipses = FALSE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Diagnosis") 
fviz_pca_var(pca.result,
             #col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
pca.ind.coord <- pca.result$x
pca.final <- cbind(PCA_, pca.ind.coord)
pca.final$kclusters <- as.character(pca.final$kclusters)

a <-fviz_pca_var(pca.result,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
)


plot (a)

# PCA Plots
# Creates plots of the first two PCs against age as well as PCA loadings vectors. These graphs would later be combined in 
#Illustrator to create nicer figures 

PC1 <- pca.final %>% ggplot(aes(PC1, Age_median)) + geom_path() + 
  geom_point(aes(PC1, Age_median, fill = kclusters), pch =21, size = 2) + 
  scale_y_reverse(breaks=seq(0,10500,500)) + 
  theme_bw() +
  scale_fill_brewer(palette = "Greens") 
PC1

PC2 <- pca.final %>% ggplot(aes(PC2, Age_median)) + geom_path() + 
  geom_point(aes(PC2, Age_median, fill = kclusters), pch =21, size = 2) + 
  scale_y_reverse(breaks=seq(0,10500,500)) + 
  theme_bw() + 
  scale_fill_brewer(palette = "Greens")
PC2

pca.final %>% ggplot(aes(PC1, PC2, fill = kclusters)) +
  geom_point(pch =21, size = 3) + 
  theme_bw() + 
  scale_fill_brewer(palette = "Greens")  

fviz_pca_var(pca.result,
             #col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
) 


