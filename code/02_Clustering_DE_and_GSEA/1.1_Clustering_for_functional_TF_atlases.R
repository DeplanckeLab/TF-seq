# title: "Clustering analysis on functional TF atlases"
# author: "Wangjie Liu"
# date: "2024/10/23"


setwd("./")

suppressPackageStartupMessages(library(Seurat)) 
suppressPackageStartupMessages(library(ggplot2)) 
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(clustree))
suppressPackageStartupMessages(library(Hmisc))


## --------------------------------------------------Figure 2b, clustering analysis on G1-phase functional TF atlas
seu <- readRDS("results/C3H10_10X_all_exps_integrated_functional_TF_atlas_G1_Phase.rds")

# find clusters for the functional TF atlas 
seu <- seu %>% FindNeighbors(dims = 1:100)
# run different resolution and then clustree to find an optimal resolution for clustering
for (i in seq(0.1, 1.5, 0.1)){
  seu <- FindClusters(seu, resolution = i, verbose = FALSE)
  p <- DimPlot(seu, group.by = "seurat_clusters",label = T)+labs(title = paste0("resolution_",i))
  print(p)
}
clustree(seu)

# apply integrated_snn_res.0.2 for Figure 2b to present general cellular processes in the atlas
cols.cluster <- c("#FFCCCC","#99CCFF","#E69F00","#CCCC00","#DDCC77","#CC99CC","#F0E442","#FF9999","#882255","#009E73","#56B4E9","darkgreen","#D55E00","#99FF99","#0072B2")
names(cols.cluster) <- as.character(0:14)
cols.cluster
p <- DimPlot(seu, group.by = "integrated_snn_res.0.2",label = T, order = T)+
  scale_color_manual(labels = c("0", "1", "2", "3", "4", "5","6","7", "8","9","10","11","12","13","14"), values = cols.cluster)+
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
p
ggsave(p, filename = "figures/Fig_2b.pdf", width = 6.5, height = 5)


## --------------------------------------------------Figure 4.1a, b; clustering analysis on G1-phase-corrected functional TF atlas

seu <- readRDS("results/C3H10_10X_all_exps_integrated_functional_TF_atlas_G1_Phase_corrected.rds")

seu <- FindNeighbors(seu, dims = 1:100)
for (i in seq(0.1, 1.5, 0.1)){
  seu <- FindClusters(seu, resolution = i, verbose = FALSE)
  p <- DimPlot(seu, group.by = "seurat_clusters",label = T)+labs(title = paste0("resolution_",i))
  print(p)
}
p.clustree <- clustree(seu, node_text_size = 0, edge_width = 1)
p.clustree
ggsave(p.clustree, filename = "figures/Fig_S4.1a.pdf",
       width = 10, height = 14, units = "in")

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
ggsave(p_umap, filename = "figures/Fig_S4.1b.pdf", width = 7, height = 5)