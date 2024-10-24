
#--------------------- CLUSTER PROFILER ---------------------

setwd("Your Enrichment Directory")
#remotes::install_github("YuLab-SMU/clusterProfiler") 
#https://www.biostars.org/p/305258/
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(DOSE)
library(cowplot)

#Load list of genes as character:
ent_gene <- degs.ag$Gene.ID

########## ClusterProfiler (GO & KEGG) ########## 
## GO enrichment 
ego <- enrichGO(gene = ent_gene,
                OrgDb = org.Hs.eg.db ,
                ont = "ALL" ,
)

#use slice function to select top 10 in each 3 GO terms:(and make a dataframe)

egom <-arrange(ego@result, p.adjust) %>% 
  group_by(ONTOLOGY) %>% 
  slice(1:10)

#cALCULATE RICH FACTOR:
egom <- mutate(egom, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

#write enrichGO
ego_table <- mutate(ego@result, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

write.table(ego_table, file="enrichGO.txt", quote=F, sep="\t")

####### dotplot BY Rich factor:

d1 <- ggplot(egom, 
             aes(richFactor,Description)) + 
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_colour_gradientn(
    colours = rainbow(10, s = 1, v = 1, start = 0, end = max(1, 10 - 1)/10,),
    name="Adjusted P-Value")+
  scale_size_continuous(range=c(2, 6)) +
  theme_bw(base_size = 14)+ 
  xlab("Rich factor") +
  ylab(NULL) +  
  guides(size=guide_legend("Gene number",order = 1),
         "Adjusted P-Value" = guide_legend(order = 2))+
  ggtitle("Statistics of GO Enrichment")+
  facet_grid(rows = vars(ONTOLOGY),scales = "free_y")


## KEGG enrichment :
options(clusterProfiler.download.method = "wininet") # cause KEGG API has moved to HTTPS 
#(https://github.com/YuLab-SMU/clusterProfiler/issues/470)

#gene_vector <- setNames(rep(1, length(ent_gene)), ent_gene)

ekg <- enrichKEGG(ent_gene,
                  organism = "hsa",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1)

#deactive pvalue and qvalue effect on the numbers of pathways in results (get all results)

view(ekg)

#view gene symbols:

ekgm<- setReadable(ekg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

#cALCULATE RICH FACTOR:
ekgm <- mutate(ekgm, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

#make data frame to read by ggplot2:
ekgm <- as.data.frame(ekgm)
#write enrichKEGG
write.table(ekgm, file="enrichKEGG.txt", quote=F, sep="\t")


####### dotplot BY Rich factor:

# select first 20:
ekgm <-arrange(ekgm, p.adjust) %>% slice(1:20)

d2 <- ggplot(ekgm,
             aes(richFactor,Description)) + 
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_colour_gradientn(
    colours = rainbow(10, s = 1, v = 1, start = 0, end = max(1, 10 - 1)/10,),
    name="Adjusted P-Value")+
  scale_size_continuous(range=c(2, 6)) +
  theme_bw(base_size = 14)+ 
  xlab("Rich factor") +
  ylab(NULL) +  guides(size=guide_legend("Gene number"))+
  ggtitle("Statistics of KEGG Enrichment")


#COMBINE DOTPLOT:
tiff("dotplot_combine.tiff",width = 5500, height = 2300,res = 300,units = "px")

plot_grid(
  d2 ,d1,
  rel_widths = c(1.75, 2),
  nrow = 1,
  labels = c('A', 'B'), 
  label_fontfamily = "serif",
  label_fontface = "bold",
  label_colour = "black",
  label_size = 30
)

dev.off()



