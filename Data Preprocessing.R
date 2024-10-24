
#---------Expression data from esophageal squamous cell carcinoma patients---------------
### Esophageal Cancer vs Paratumor (CGC Congress 2023- www.cgc2023.com) ###
#GSE161533  #Samples (84 total) 
#28 Paratumor and 28 Tumor specimens
#GPL570	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array

#Load Packages:
# Load Required Packages
library(Biobase)       # Bioinformatics data structure
library(limma)         # Linear Models for Microarray Data
library(GEOquery)      # Interface to the NCBI GEO data
library(pheatmap)      # For creating heatmaps
library(gplots)        # Additional plotting functions
library(ggplot2)       # For advanced visualizations
library(reshape2)      # Reshaping data for visualization
library(plyr)          # Data manipulation

#Load previous data to speed up:
setwd("Your Directory")
lapply(X = c("DEGs.RData","EXDATA.RData","HEATMAPOO.RData"),load,.GlobalEnv)

#--------------------- LOAD DATA -------------------------
setwd("Your Directory")

series <- "GSE39144" 
gset <- getGEO(series, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "dataq/")
platform<- "GPL570"
if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
#Export Expression Matrix (Expression data) :
ex <- exprs(gset)
dim(ex)


#--------------------- SAMPLE SELECTION ---------------------
#Sample Selection:
gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXX1111111111111111111111",
               "1111110000000000000000000000000000")
sml <- strsplit(gsms, split="")[[1]]

# X = Normal  1 = paratumor 0 = tumor

# filter out excluded samples (marked as "X")
sel <- which(sml != "X") #here you can choose which sample you want to exclude
sml <- sml[sel]
gset <- gset[ ,sel]
#Rewrite expression for selected samples:
##in this stage we define ex without norma samples
ex <- exprs(gset)
dim(ex)

#Log2 transformation needed?
max(ex)
min(ex)
#Yes or No???
## if YES:
ex<- log2(ex+1) #The data is scaled!
# Apply ITTTTT!!!
exprs(gset) <- ex

#for plots:
gr <- factor(c(rep("Paratumor",28), rep("Tumor", 28)))

#Save Data:
setwd("Your Directory")
save(data, ex, gr, file = "EXDATA.RData")

#--------------------- QC PLOTS ---------------------
setwd("E:\\Projects\\CGC2023/CGC Mammad_Esophageal Cancer/plots/")

### Draw BoxPlot:
pdf("boxplot.pdf")
boxplot(ex)
dev.off()


### Draw PCA:
#Method 1
exm <- ex-rowMeans(ex) #rowMeans is a vector 
pc <- prcomp(exm)
pcr <-data.frame(pc$r[,1:3], gr)


png("PCA_rowMeans_N vs C.png",width = 2000, height = 1500,res = 300,units = "px")

ggplot(pcr,aes(PC1,PC2,color=gr)) + ggtitle("PCA rowMeans")+
  geom_label(size=2.2, fontface = "bold",label= gr)+
  geom_vline(xintercept=0, col="black",linetype = "longdash") + 
  geom_hline(yintercept=0, col="black",linetype = "longdash") + 
  stat_ellipse(aes(fill=gr),type = "norm",geom = "polygon",
               alpha = 0.09,level = 0.6)+
  theme_bw()

dev.off()

#Method 2
ex.scale <- t(scale(t((ex)), scale=T,  center = F))
pc <- prcomp(ex.scale)
pcr <- data.frame(pc$r[,1:3], Group=gr)


png("PCA_scale_N vs C.png",width = 2000, height = 1500,res = 300,units = "px")

ggplot(pcr,aes(PC1,PC2,color=gr)) + ggtitle("PCA scaled")+
  geom_label(size=2.2, fontface = "bold",label= gr)+
  geom_vline(xintercept=0, col="black",linetype = "longdash") + 
  geom_hline(yintercept=0, col="black",linetype = "longdash") + 
  stat_ellipse(aes(fill=gr),type = "norm",geom = "polygon",
               alpha = 0.09,level = 0.6)+
  theme_bw()


dev.off()

