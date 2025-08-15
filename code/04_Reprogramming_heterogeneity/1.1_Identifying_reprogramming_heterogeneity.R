# title: "Identifying TFs driving reprogramming heterogeneity"
# input data: "functional TF atlas"
# author: "Wangjie Liu"

setwd("./")


suppressPackageStartupMessages(library(Seurat)) 
suppressPackageStartupMessages(library(ggplot2)) 
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(clustree)) 

## read functional TF atlas with functional TF cells, D0 and ref cells in G1 at adjusted phase
seu <- readRDS("results/C3H10_10X_all_exps_integrated_functional_TF_atlas_G1_Phase_corrected.rds")

## -------------------------------------------------- Extended_data_figure_5a
set.seed(10)
# find clusters for the functional TF atlas
# run different resolution and then clustree to find an optimal resolution for clustering
seu <- FindNeighbors(seu, dims = 1:100)
for (i in seq(0.1, 1.5, 0.1)){
  seu <- FindClusters(seu, resolution = i, verbose = FALSE)
  p <- DimPlot(seu, group.by = "seurat_clusters",label = T)+labs(title = paste0("resolution_",i))
  print(p)
}

p.clustree <- clustree(seu, node_text_size = 0, edge_width = 1)
p.clustree
ggsave(p.clustree, filename = "figures/EDFig_5a.pdf",
       width = 10, height = 14, units = "in")



## -------------------------------------------------- Extended_data_figure_5b 
# choose resolution 1.2 for clustering, because the number of clusters becomes stable 
seu$seurat_clusters <- seu$integrated_snn_res.1.2
p_umap <- DimPlot(seu, group.by = "seurat_clusters", order = T, label = F, cols = Cols_res_0.6)+
  theme(axis.line = element_blank(),
        legend.position = "none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.background=element_blank(),
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = NULL)
p_umap
ggsave(p_umap, filename = "figures/EDFig_5b.pdf", width = 7, height = 5)


## -------------------------------------------------- Extended_data_figure_5c
### extract clustering metadata at "integrated_snn_res.1.2" and find TFs with more than 1 cluster
df <- data.frame(row.names = rownames(seu@meta.data),
                 TF = seu$TF,
                 Cluster = seu$integrated_snn_res.1.2)
df$Cluster <- as.factor(df$Cluster)
df <- df %>% group_by(TF) %>% mutate(Number_of_TF_cells = length(TF))
df <- df[df$Number_of_TF_cells >= 30,] 
df <- lapply(unique(df$Cluster), function(x){
  df <- df[df$Cluster == x, ] %>% group_by(TF) %>% mutate(Number_of_TF_cells_in_each_clusters = length(TF))
  df <- df[df$Number_of_TF_cells_in_each_clusters >= 10,] # minimum TF cells in each cluster
  return(df)
})
df <- data.table::rbindlist(df)
df_percent <- lapply(unique(df$Cluster), function(x){
  Percentage_TF <- df[df$Cluster == x, ] %>% group_by(TF) %>% mutate(percentage_TF = length(TF)/Number_of_TF_cells)
  Percentage_TF <- Percentage_TF[!duplicated(Percentage_TF),]
  return(Percentage_TF)
})
df_percent <- data.table::rbindlist(df_percent)
df_percent$percentage_TF <- df_percent$percentage_TF*100 %>% as.numeric()
df_percent$percentage_TF <- round(df_percent$percentage_TF)
df_percent <- df_percent[df_percent$percentage_TF > 5 & df_percent$percentage_TF <= 95,]
D0.clusters <- df_percent$Cluster[df_percent$TF == "D0" & df_percent$percentage_TF >= 15] %>% as.character() # main clusters of D0 (non-confluent only)
D0.confluent.clusters <- df_percent$Cluster[df_percent$TF == "D0_confluent" & df_percent$percentage_TF >= 15] %>% as.character()
TF.keep <- table(df_percent$TF[!df_percent$Cluster %in% D0.confluent.clusters]) 
TF.keep <- TF.keep[TF.keep > 1] %>% names()
df_percent <- df_percent[df_percent$TF %in% c(TF.keep, "D0","D0_confluent"),]

nCluster <- unique(df_percent$Cluster) %>% as.character()
nCluster <- nCluster[order(as.numeric(nCluster), decreasing = F)]

# reform the data frame for heatmap plot
df_mat <- lapply(unique(df_percent$TF), function(x){
  df.reform <- data.frame(row.names = df_percent$Cluster[df_percent$TF == x],
                          percentage_TF = df_percent$percentage_TF[df_percent$TF == x])
  Cluster.supplement <- nCluster[!nCluster %in% rownames(df.reform)]
  df.reform <- rbind(df.reform, data.frame(row.names = Cluster.supplement,
                                           percentage_TF = rep(0, length(Cluster.supplement)))) 
  df.reform <- t(df.reform) %>% as.data.frame()
  df.reform$TF <- x
  df.reform <- df.reform[, c(as.character(nCluster),"TF")]
  return(df.reform)
})

df_mat <- data.table::rbindlist(df_mat) 
rownames.df_mat <- df_mat$TF
df_mat$TF <- NULL
df_mat <- as.matrix(df_mat)
rownames(df_mat) <- rownames.df_mat

#remove TF combination, adipo and myo ref cells
TF.com <- c("Cebpa-Mycn","Cebpa-Myog","Cebpa-Pparg","Mycn-Myog","Mycn-Pparg","Mycn-Runx2","Pparg-Runx2")
df_mat <- df_mat[!rownames(df_mat) %in% c(TF.com, "Adipo_ref","Myo_ref"),]

rownames(df_mat) <- toupper(rownames(df_mat))
rownames(df_mat)[rownames(df_mat) == "D0"] <- "Ctr.non.conf"
rownames(df_mat)[rownames(df_mat) == "D0_CONFLUENT"] <- "Ctr.conf"
p_heatmap <- pheatmap::pheatmap(df_mat, 
                                cellheight = 10,
                                cellwidth = 12, 
                                fontsize_row = 10, 
                                fontsize_col = 10,
                                color = colorRampPalette(c("#edf8e9","#006d2c"))(100),
                                border_color = NA,
                                cluster_rows = T, 
                                clustering_method = "ward.D2",
                                cluster_cols = F,
                                treeheight_row = 0)

p_heatmap

ggsave(p_heatmap, filename = "figures/EDFig_5c.pdf",
       width = 8, height = 12)