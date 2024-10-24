
#--------------------- Find DEGs ---------------------
setwd("Your Directory")
getwd()
## Make Design Matrix:
# assign samples to groups and set up design matrix
gs <- factor(sml) #   0 = Cancer , 1 = Normal
groups <- make.names(c("Tumor","Paratumor")) #the order of the names should be by Levels(gs) order
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

#Limma starts:
fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-") 
cts #check: "Tumor-Paratumor" / "Cancer-Normal" / "Experiment-Control"

cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust.method ="fdr", sort.by="logFC", number=Inf)
tT <- subset(tT, select=c("Gene.symbol","Gene.ID","logFC","adj.P.Val","P.Value","t","B","ID"))

head(tT)

#remove non gene probes:
tT <- tT[-which(tT$Gene.symbol == ""), ]

#remove /// in gene symbols and select first name:
tT$Gene.symbol <- sub("///.*","",tT$Gene.symbol)
#remove /// in gene IDs and select first ID:
tT$Gene.ID <- sub("///.*","",tT$Gene.ID)

#write.table(tT, file="total results.txt", quote=F, sep="\t")

# Total DEGs(up and down): for main heatmap

#DEGs:#  IDs
degs <- subset(tT, abs(logFC) > 2 & adj.P.Val < 0.05)

###  multiple probe sets mapping to the same gene: A need for aggregation (Thanks to Dr. Bahari) ###

nums <- table(degs[,"Gene.symbol"])
degs.ag <- data.frame(matrix(0,nrow = length(nums), ncol = ncol(degs))) 

for(i in 1:length(nums)){
  j <- which(degs[,"Gene.symbol"] == names(nums[i]))
  fac <- degs[j,"logFC"] #(cause we want max LogFC)
  k <- which(max(fac) == fac) #"max" can be changed by the other functions (e.g. mean, min, median and etc.)
  degs.ag[i,] <- degs[j[k],]
}

colnames(degs.ag) <- colnames(degs)
rownames(degs.ag) <- degs.ag[,"Gene.symbol"]

#Find how much upregulate and downregulate genes you had:
degsup <- subset(degs.ag, logFC > 2 & adj.P.Val < 0.05) #OncoGene
degsdown <- subset(degs.ag, logFC < -2 & adj.P.Val < 0.05) #TumorSUpressor

#"degs.ag" is what we need for further analysis:
write.table(degs.ag, file="Degs - aggregated in R.txt", quote=F, sep="\t",
            row.names = F)

#Write down gene symbols for STRING input:
write.table(degs.ag$Gene.symbol, file="Degs GeneNames.txt", quote=F,
            row.names = F, col.names = F)

#Save Data:
setwd("Plots Directory")
save(degs,degs.ag,degsdown,degsup,design,tT, file = "DEGs.RData")
