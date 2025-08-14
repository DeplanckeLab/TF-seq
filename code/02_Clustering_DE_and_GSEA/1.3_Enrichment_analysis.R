# title: "Pathway or Gene set enrichment analysis on G1-phase functional TF atlas"
# author: "Wangjie Liu"
# date: "2024/10/23"


setwd("./")




library(Seurat)
library(tidyverse)
library(dplyr)
library(gridExtra)
library(fgsea)
library(enrichplot)
library(ggplot2)
library(tibble) 


# load gene annotation
data.annot <- read.table("data/GRCm38.96_Vector_data.annot.txt", sep = "\t")
rownames(data.annot) <- data.annot$ens_id

# load the seurat object of G1-phase functional TF atlas
seu <- readRDS("C3H10_10X_all_exps_integrated_functional_TF_atlas_G1_Phase.rds")

# filter out genes that are not detected in at least 3 cells, to not diminish GSEA signals
features.filtered <- rowSums(seu@assays$RNA@counts > 0) >=10
features.filtered <- names(features.filtered)[features.filtered == TRUE]

# prepare Panglao gene markers for cell types of interest
Panglao <- readRDS("data/PanglaoDB_markers_REFORMAT_15_Sep_2019.Rds")
f <- function(x){
  x <- x[x$Mm,] # select mouse genes (Mm = TRUE)
  x <- x[order(x$specificity_mouse, decreasing = T), "ens_id"]
  x <- x[!is.na(x)]
  return(x)
}
Panglao_perCelltypes <- lapply(Panglao, f)

# select lineages of interest (MSC)
lineage.ois <- list(Panglao_perCelltypes$Adipocytes, Panglao_perCelltypes$Myoblasts, Panglao_perCelltypes$Chondrocytes, Panglao_perCelltypes$Osteoblasts)
names(lineage.ois) <- c("Adipocytes","Myoblasts", "Chondrocytes","Osteoblasts") 

# filter the Panglao gene markers
lineage.ois <- lapply(lineage.ois, function(x){
  x <- x[x %in% features.filtered]
  return(x)
})

# customized adipo markers
adipo_markers <- c("Fabp4", "Lpl", "Pparg", "Lipe", "Adipoq", "Cd36", "Plin4", "Plin2", "Plin1", "Cebpa", "Cebpb", "Cidec", "Cidea")
adipo_markers_id <- data.annot[data.annot$gene_short_name %in% adipo_markers,"ens_id"]

# convert between ensembl id or gene symbol as needed
ens_id <- F
if (ens_id == F){
  ## for gene symbol
  #  convert ensembl id to gene symbol
  lineage.ois <- lapply(lineage.ois, function(x){x <- data.annot$gene_short_name[data.annot$ens_id %in% x]})
  # replace adipo lineage of Panglao by customized markers representing more mature and specific adipose markers
  lineage.ois$Adipocytes <-  adipo_markers
  
} else {
  ## for ens_id
  #  convert gene symbol of customized adipo markers to ens_id
  lineage.ois$Adipocytes <- adipo_markers_id
}

# Load DEGs of all TFs
DEGs.TFs <- read.csv("results/DEGs_glmQLFtest_by_TF_batchoverall_allD0_functionalCells.rds")
DEGs.TFs <- DEGs.TFs[DEGs.TFs$ens_id %in% features.filtered,]

# TFs that have less than 25 cells in total are filtered out 
nCell.tfois <- table(seu$TF)
TFois <- nCell.tfois[nCell.tfois >=25 ] %>% names() 

# Exclude controls, TF combinations
TFois <- TFois[!TFois %in% c("Pparg-Runx2",
                             "Mycn-Runx2",
                             "Mycn-Myog",
                             "Mycn-Pparg",
                             "Cebpa-Mycn",
                             "Cebpa-Myog",
                             "Cebpa-Pparg",
                             "D0",
                             "D0_confluent")] 


## ------------------------------------------- run fgsea on Panglao genesets
DEGs.TFs <- DEGs.TFs[DEGs.TFs$TF %in% TFois, ] 
table(DEGs.TFs$TF) %>% length() # 101 TFs
Padj_threshold <- 0.1

df.gsea <- lapply(unique(DEGs.TFs$TF), function(x){
  DEG.oi <- DEGs.TFs[DEGs.TFs$TF == x,]
  DEG.oi <- DEG.oi[order(DEG.oi$logFC, decreasing = T), ]
  ranks <- DEG.oi$logFC
  names(ranks) <-  DEG.oi$ens_id
  set.seed(10)
  fgseaRes <- fgsea(pathways = lineage.ois, stats = ranks)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>% filter(padj < Padj_threshold)

  fgseaResTidy$TF <- x
  if(nrow(fgseaResTidy) == 0) {
    fgseaResTidy <- NA
  }
  
  return(fgseaResTidy)
})


df.gsea <- df.gsea[!is.na(df.gsea)]
df.gsea <- data.table::rbindlist(df.gsea, use.names=TRUE) %>% as.data.frame()
df.gsea$TF <- as.factor(df.gsea$TF)
df.gsea <- df.gsea[order(df.gsea$TF, decreasing = T),]
df.gsea$`-log10(Padj)` <- -log10(df.gsea$padj)


## ------------------------------------------- plotting gsea results
# clustering/ordering TFs by pathways
df.gsea$pathway <- as.factor(df.gsea$pathway)
TFs.gsea.all <- unique(df.gsea$TF) %>% as.vector()

# TFs that are positively enriched for at least one pathway or TFs that negatvely enriched for all pathways
TFs.gsea.Pos <- unique(df.gsea$TF[df.gsea$NES >0]) %>% as.vector()
TFs.gsea.Neg <- TFs.gsea.all[!TFs.gsea.all %in% TFs.gsea.Pos]

df.order<- lapply(TFs.gsea.Pos, function(x){
  df.gsea.sub <- df.gsea[df.gsea$TF == x,]
  main.pathway <- df.gsea.sub$pathway[which.max(df.gsea.sub$NES)] %>% as.vector()
  df.mainES <- data.frame(TF = x, main_ES = main.pathway)
  return(df.mainES)
})
df.order <- data.table::rbindlist(df.order) %>% as.data.frame()
df.order <- df.order[order(df.order$main_ES, decreasing = F),]
TF_order <- df.order$TF

df.gsea$pathway <- factor(df.gsea$pathway, levels = c("Adipocytes","Chondrocytes","Myoblasts","Osteoblasts"))
df.gsea$TF <- factor(df.gsea$TF, levels = c(TF_order, TFs.gsea.Neg))


## for horizontal; 
p <- ggplot(df.gsea)+geom_point(aes(x = TF, y = pathway, size = padj, col = NES))+
  viridis::scale_color_viridis()+
  scale_size(name = expression(italic("p.adj")), trans = "reverse")+
  mashaGgplot2Theme+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, color = "black"),
        axis.text.y = element_text(face = "bold", color = "black"),
        axis.title.y = element_blank())
p


ggsave(p, filename = "figures/Figure2d",
       width = 16, height = 4.5)  


##------------------------prepare ranked genelist for GSEA software (which generated Extended_data_figure_3d with the input of ranked genelist)

## rank genes based on DE
DE.cluster.oi <- readRDS("results/DEGs_glmQLFtest_by_Cluster_8_VS_Clusters_Ctr_at_R0.2_G1_Phase_batchcovariate.rds")
DE.cluster.oi <- DE.cluster.oi[order(DE.cluster.oi$logFC, decreasing = T), ]
rank.oi <- DE.cluster.oi$logFC
names(rank.oi) <-  DE.cluster.oi$gene_short_name
geneList.ranked <- data.frame(names(rank.oi), as.numeric(rank.oi))
write.table(geneList.ranked, file = "results/Preranked_genelist.txt", col.names = F, row.names = F, sep = "\t", quote = F)



