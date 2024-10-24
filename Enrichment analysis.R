#--------------------- Downstream Analysis ---------------------
#Leze bOoOo:


#--------------------- VolcanoPlots for DEGs ---------------------

setwd("Plots Directory")
volc <- subset(tT)

volc$Significant <- "No"
volc$Significant[volc$logFC > 2 & volc$adj.P.Val < 0.05] <- "Up"
volc$Significant[volc$logFC < -2 & volc$adj.P.Val < 0.05] <- "Down"

#tiff("volcano plot.tiff",height = 1600, width = 2000,res=300,units = "px" )
pdf("volcano.pdf")

ggplot(volc, aes(logFC, -log10(adj.P.Val),color =Significant))+ 
  geom_point(size = 1.6,shape = 19) + theme_bw()+ 
  geom_vline(xintercept=c(-2, 2), col="#EF00FF",linetype = "longdash") + 
  geom_hline(yintercept=-log10(0.05), col="#EF00FF",linetype = "longdash") +
  scale_color_manual(values=c("blue", "gray", "red")) +
  ggtitle("Volcano Plot")



ggplot(volc, aes(logFC, -log10(adj.P.Val), color = Significant)) + 
  geom_point(size = 1.6, shape = 19) + 
  theme_bw() + 
  geom_vline(xintercept=c(-2, 2), col="#EF00FF",linetype = "longdash") + 
  geom_hline(yintercept=-log10(0.05), col="#EF00FF",linetype = "longdash") +
  scale_color_manual(values=c("blue", "gray", "red")) +
  ggtitle("Volcano Plot") +
  geom_text(aes(label = Gene.symbol), size = 1.5, nudge_y = 0.2)


dev.off()


#--------------------- HEATMAP for DEGs ---------------------
#-------------- Heatmapoo is an expression matrix of my DEGs (also for WGCNA)
#--- NEW WAY! ---#
#make a dataframe of expression matrix:
bayan <- data.frame(ex)
# extract expression of DEGs for downstream analysis
heatmapoo <- bayan[rownames(bayan) %in% degs.ag$ID,]

#replace probe IDs in rownames heatmapoo into gene symbol:
# create a named vector of gene symbols from degs.ag
symbol_match <- setNames(degs.ag$Gene.symbol, degs.ag$ID)
# replace rownames in heatmapoo with gene symbols
rownames(heatmapoo) <- symbol_match[rownames(heatmapoo)]

#Save Data:
setwd("Your Directory")
save(heatmapoo, file = "HEATMAPOO.RData")
#--------------#

#making annotation_col:
coloo <- colnames(heatmapoo) 
annotation1 <- data.frame(gr, coloo)
rownames(annotation1) <- annotation1$coloo
colnames(annotation1) = "Group"
annotation1 <- annotation1[,-2,drop=F]


annotation2 <- data.frame(degs.ag$logFC)
rownames(annotation2) <- degs.ag$Gene.symbol
colnames(annotation2) = "Up/Down regulated"
annotation2$`Up/Down regulated`[degs.ag$logFC > 2 & degs.ag$adj.P.Val < 0.05] <- "Up"
annotation2$`Up/Down regulated`[degs.ag$logFC < -2 & degs.ag$adj.P.Val < 0.05] <- "Down"


ann_colors = list(
  "Group" = c( Paratumor= "#7570B3",  Tumor= "#3CF20F"),
  "Up/Down regulated" = c(Down = "#D7A3E6", Up = "#E7298A"))
#get the color code from rapidtables RGB colors
#https://www.rapidtables.com/web/color/RGB_Color.html

setwd("Your Plot Directory")

png("heatmap_default_2.png",height = 1600, width = 2000,res=300,units = "px" )

pheatmap(heatmapoo,color= bluered(256), show_rownames=F, show_colnames = F,
         border_color =NA,annotation_colors = ann_colors,
         annotation_col = annotation1, annotation_row = annotation2, 
         scale = "row" , angle_col=45,fontsize_col=4.5, fontsize_row = 4.5,
         annotation_names_row = F,annotation_names_col = F,
         clustering_distance_rows = "correlation", #Pearson correlation
         clustering_distance_cols= "correlation",  #Pearson correlation
         clustering_method = "complete")


dev.off()

